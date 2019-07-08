//This file is part of faust, faster annotation using shape-constrained trees.     

//faust is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.

//faust is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.

//You should have received a copy of the GNU General Public License
//along with faust.  If not, see <http://www.gnu.org/licenses/>.
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <Rcpp.h>
#include "faust.h"
#include <stack>
#include <random>
#include <map>
#include <set>


//The search for candidate clusters occurs over the space of all clustering trees.
//The search begings by selecting a multimodal column from the input dataMatrix.
//We represent dataMatrix as std::vector<std::vector<double>> to emphazies that the method views the dataMatrix
//as a related colection of 1-d column vectors of the same length; the method does not compute any matrix quantities (SVD, inverse, PCA, etc.)
//The selection occurs using singleDip, which produces a p-value for each column of the input dataMatrix.
//Any column vector selected consistutes the root of a clustering tree.
//Modal subsets of the selected column are induced using tsGates with modePrior=0(primarily).
//(for small subsets, doubleDip and tripleDip are used to set modePrior to non-zero values.)
//The modal subsets are then collected onto a stack of subsetInfo, which we define here.
struct subsetInfo {
  std::vector<bool> activeRowIndices; //true if the observation is in the current subset
  std::vector<bool> activeColumnIndices; //true if the column of the data matrix is viable for recursion.
  int nodeDepth;
  double nodePathScore;
};
//Each subset on the stack is then examined for multimodality in turn.
//A candidate cluster is any subset of the dataMatrix which is unimodal along all columns.

/*
unsigned long long sizeOfGateInfo(gateInfo& gateInfoRef) {
  unsigned long long totalSize = 0;
  totalSize += sizeof(int);
  totalSize += sizeof(unsigned long);
  std::vector<double> agates = gateInfoRef.gates;
  totalSize += (agates.size() * sizeof(double));
  return totalSize;
}
*/

struct unrestrictedSubsetData {
  std::vector<std::vector<double>> urSubset;
  std::vector<long> urCountVec; //count how many observations are in each urSubset vector
};
  
unrestrictedSubsetData getActiveDataSubset(const std::vector<bool>& activeRowsInSubset,
					   const std::vector<bool>& activeColsInSubset,
					   const std::vector<std::vector<double>>& fullDataMatrix)
{
  //obtain active subset of data matrix in currentDataMatrix.
  std::vector<std::vector<double>> activeMatrix;
  std::vector<double> tmpC, tmpR;
  std::vector<long> unrestrictedCountVector;
  long tmpCount;
  for (int i =0; i != activeColsInSubset.size(); i++) {
    tmpC.clear();
    tmpR.clear();
    tmpC = fullDataMatrix[i];
    tmpCount = 0;
    for (auto j = 0; j != activeRowsInSubset.size(); j++) {
      if (activeRowsInSubset[j] == true) {
	tmpR.push_back(tmpC[j]);
	++tmpCount;
      }
    }
    activeMatrix.push_back(tmpR);
    unrestrictedCountVector.push_back(tmpCount);
  }
  unrestrictedSubsetData retData = {activeMatrix, unrestrictedCountVector};
  return retData;
}

//if the search is restricted, we use the following structure when working with unrestricted subsets
struct allSubsetData {
  std::vector<std::vector<double>> fullSubset;
  std::vector<std::vector<double>> urSubset;
  std::vector<bool> lowResV;
  std::vector<bool> highResV;
  std::vector<long> urCountVec; //count how many observations are in each urSubset vector
  std::vector<std::vector<long>> restrictedIndices;
};

allSubsetData getAllDataSubsets(const std::vector<bool>& activeRowsInSubset,
				const std::vector<bool>& activeColsInSubset,
				const std::vector<std::vector<double>>& fullDataMatrix,
				const std::vector<std::vector<int>>& resMat,
				bool recordCounts)
{
  std::vector<std::vector<double>> fullDataSubset, unrestrictedDataSubset;
  std::vector<double> tmpC, tmpR,tmpU;
  std::vector<int> tmpRestrict;
  //the following vectors indicate if a column *ever* has a restriction of the stated type.
  std::vector<bool> lowRestrictionVector, highRestrictionVector;
  bool lowRes, highRes;
  std::vector<long> unrestrictedCountVector;
  std::vector<std::vector<long>> restrictedIndexMat;
  std::vector<long> restrictedIndexVector;
  long tmpCount;
  //obtain active subset of data matrix in currentDataMatrix.
  //record the cooresponding rows in the 
  for (int i =0; i != activeColsInSubset.size(); ++i) {
    tmpC.clear();
    tmpR.clear();
    tmpU.clear();
    tmpRestrict.clear();
    restrictedIndexVector.clear();
    restrictedIndexVector.clear();
    tmpC = fullDataMatrix[i];
    tmpRestrict = resMat[i];
    lowRes = false;
    highRes = false;
    tmpCount = 0;
    for (auto j = 0; j != activeRowsInSubset.size(); ++j) {
      if (activeRowsInSubset[j] == true) {
	tmpR.push_back(tmpC[j]); //tmpR contains all observations in the active set.
	if (tmpRestrict[j] == 0) {
	  tmpU.push_back(tmpC[j]); //tmpU contains only unrestricted observations in the active set.
	  ++tmpCount;
	  if (recordCounts) {
	    restrictedIndexVector.push_back(j);
	  }
	}
	else if (tmpRestrict[j] == 1) { //user signal for low restriction
	  lowRes = true;
	}
	else if (tmpRestrict[j] == 2) { //user signal for high restriction.
	  highRes = true;
	}
      }
    }
    fullDataSubset.push_back(tmpR);
    unrestrictedDataSubset.push_back(tmpU);
    lowRestrictionVector.push_back(lowRes);
    highRestrictionVector.push_back(highRes);
    unrestrictedCountVector.push_back(tmpCount);
    if (recordCounts) {
      restrictedIndexMat.push_back(restrictedIndexVector);
    }
  }

  if (!recordCounts) {
    restrictedIndexVector.push_back(0);
    restrictedIndexMat.push_back(restrictedIndexVector);
  }
  
  allSubsetData subDat = {fullDataSubset, unrestrictedDataSubset, lowRestrictionVector,
			  highRestrictionVector, unrestrictedCountVector,restrictedIndexMat};

  return subDat;
}


struct viableColumnData {
  std::vector<int> columnIndices;
  std::vector<double> columnScores;
};

viableColumnData getViableColumns(const std::vector<std::vector<double>>& sortedData,
				  const std::vector<bool>& actCols,
				  double singleDipThreshold,
				  bool randomSearch,
				  std::mt19937& twRef,
				  bool parexSearch,
				  bool& firstParexIter,
				  int parexColNum)
{
  std::vector<int> viableColumns, returnColumns;
  std::vector<double> viableColumnScores, returnScores;
  double colDipScore;
  //along each active column, check if DIP test p-value below threshold.
  //std::cout << "Depth Scores at construction." << std::endl;
  for (int i =0; i != actCols.size(); ++i) {
    if (actCols[i] == true) {
      colDipScore = singleDip(sortedData[i]);
      if (colDipScore < singleDipThreshold) {
	viableColumns.push_back(i);
	viableColumnScores.push_back(colDipScore);
	//std::cout << "col " << i << ": " << colDipScore << std::endl;
      }
    }
  }
  //if randomSearch = true, user wants to sample from the space of clustering trees
  //pick a column uninformly at random from the viable columns, and return it.
  if ((randomSearch) && ((viableColumns.size()) > 0)) {
    auto numViable = viableColumns.size();
    //std::random_device randomDevice;  //used to obtain a seed for the random number engine
    //std::mt19937 generate(randomDevice()); //standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> columnDistribution(0, (numViable-1)); 
    int columnChoice = columnDistribution(twRef); //we select the column here
    returnColumns.push_back(viableColumns[columnChoice]);
    returnScores.push_back(viableColumnScores[columnChoice]);
  }
  else if ((parexSearch) && (firstParexIter) && ((viableColumns.size()) > 0))  {
    firstParexIter = false;
    returnColumns.push_back(parexColNum);
    colDipScore = singleDip(sortedData[parexColNum]);
    //std::cout << "parcol " << parexColNum << ": " << colDipScore << std::endl;
    returnScores.push_back(colDipScore);
  }
  else {
    //if the search is exhaustive, or parallel exhaustive and past first iter, return all viable columns.
    returnColumns = viableColumns;
    returnScores = viableColumnScores;
  }
  viableColumnData vcInfo = {returnColumns, returnScores};
  return vcInfo;
}

std::vector<double> getCutPointsForViableColumn(int currentColumn,
						const std::vector<std::vector<double>>& sdMat,
						double singleDipThreshold){
  int modeEst = 0;
  std::vector<double> sortedData = sdMat[currentColumn];
  double ddPval, tdPval;
  //If we are considering smaller sets of observations, use the double & triple dip.
  if (sortedData.size() <= 500) {
    ddPval = doubleDip(sortedData);
    if (ddPval >= (singleDipThreshold/2.0)) {
      modeEst = 2;
    }
    else if (sortedData.size() <= 400) {
      tdPval = tripleDip(sortedData);
      if (tdPval >= (singleDipThreshold/3.0)) {
	modeEst = 3;
      }
      else {
	modeEst = 4;
      }
    }
    else {
      //here we have somewhere between 400 and 500 observations.
      //since we currently lack critical values for the tripleDip in this subset size, we assume three modes
      modeEst = 3;
    }
  }
  std::vector<double> cutPoints;
  cutPoints = tsGates(sortedData,modeEst);
  return cutPoints;
}





std::vector<std::vector<double>> sortDataSubset(const std::vector<std::vector<double>>& unsortedDataSubset) {
  //the dip and taut-string require the column vectors in the active data subset to be sorted for their
  //computation. Each iteration of the candidate cluster search makes a sorted copy of the active data subset once
  //and passes the sorted subset to the dip and taut-string functions.
  std::vector<std::vector<double>> dataSubset = unsortedDataSubset;
  for (auto it = dataSubset.begin(); it != dataSubset.end(); it++){
    std::sort((*it).begin(),(*it).end());
  }
  return dataSubset;
}

void updateStackWithViableCols(const std::vector<int>& vCols,
			       const std::vector<bool>& activeColumnsInSubset,
			       const std::vector<bool>& activeRowsInSubset,
			       const std::vector<std::vector<double>>& currentDataMat,
			       const std::vector<std::vector<double>>& sortedDataMat,
			       const std::vector<long>& baseIndex,
			       std::vector<gateInfo>& gdvRef,
			       std::stack<subsetInfo>& sstackRef,
			       double dThreshold,
			       unsigned long long maximumGateNum,
			       int nodeDepth,
			       bool repeatedSplittingAllowed,
			       const std::vector<double>& vColScores,
			       double pPathScore,
			       const std::vector<long>& nodePopCounts,
			       bool recordIndices) 
{
  struct subsetInfo pushNode;
  std::vector<double> currentDataVector,cutPoints;
  std::vector<bool> rRows, rCols;
  double gateLow, gateHigh;
  gateInfo currentGateInfo;
  double childPathScore;
  std::vector<bool> dummyBool;
  //at least one column is viable. continue search for each viable column.
  for (auto i = 0; i != vCols.size(); ++i) {
    currentDataVector = currentDataMat[vCols[i]];
    //std::cout << "Getting cut points." << std::endl;
    cutPoints = getCutPointsForViableColumn(vCols[i],sortedDataMat,dThreshold);
    rCols = activeColumnsInSubset;
    if (repeatedSplittingAllowed) {
      //the user has indicated that a candidate cluster search can search a long a column vector multiple times.
      for (auto v = 0; v != rCols.size(); ++v) {
	rCols[v] = true;
      }
    }
    rCols[vCols[i]]=false; //prevent searching along same column twice.
    if (cutPoints.size() > maximumGateNum)
      continue;
    else {
      //std::cout << "Updating stack." << std::endl;
      //record information about intermediate gate locations for later processing.
      if (recordIndices) {
	//this is to record indices for later conditional plots in R.
	dummyBool.clear();
	for (auto b : activeRowsInSubset)
	  dummyBool.push_back(b);
      }
      else {
	dummyBool.clear();
	dummyBool.push_back(false);
      }
      currentGateInfo = {vCols[i], cutPoints.size(), cutPoints, nodeDepth, pPathScore,nodePopCounts[vCols[i]],dummyBool};
      gdvRef.push_back(currentGateInfo);
      childPathScore = (pPathScore * (1-vColScores[i]));
      //finally, update stack with new subsets.
      for (auto j = 1; j != cutPoints.size(); ++j) {
	gateLow = cutPoints[(j-1)];
	gateHigh = cutPoints[j];
	rRows = activeRowsInSubset;
	//find data that falls between gates. note here we read from the unsorted active data, so that
	//we recover the observation's index in the underlying data structure.
	for (auto k = 0; k != currentDataVector.size(); ++k) {
	  if ((currentDataVector[k] <= gateLow) || (currentDataVector[k] > gateHigh)) {
	    rRows[baseIndex[k]] = false;
	  }
	  if ((j == 1) && (currentDataVector[k] == gateLow)) {
	    rRows[baseIndex[k]] = true; //the activeRow minimum of current activeCol gets assigned to the first sub-collection.
	  }
	}
	pushNode = {rRows, rCols, nodeDepth, childPathScore};
	sstackRef.push(pushNode);
      }
    }
  }
  return;
}

void restrictedStackUpdate(const std::vector<int>& vCols,
			   const std::vector<bool>& activeColumnsInSubset,
			   const std::vector<bool>& activeRowsInSubset,
			   const std::vector<std::vector<double>>& currentDataMat,
			   const std::vector<std::vector<double>>& sortedDataMat,
			   const std::vector<long>& baseIndex,
			   std::vector<gateInfo>& gdvRef,
			   std::stack<subsetInfo>& sstackRef,
			   double dThreshold,
			   unsigned long long maximumGateNum,
			   int nodeDepth,
			   const std::vector<bool>& lowRestrictionVector,
			   const std::vector<bool>& highRestrictionVector,
			   bool repeatedSplittingAllowed,
			   const std::vector<double>& vColScores,
			   double pPathScore,
			   const std::vector<long>& nodePopCounts,
			   bool recordIndices) 
{
  struct subsetInfo pushNode;
  std::vector<double> currentDataVector,cutPoints;
  std::vector<bool> rRows, rCols;
  double gateLow, gateHigh;
  gateInfo currentGateInfo;
  bool columnHasLowRes, columnHasHighRes;
  double childPathScore;
  std::vector<bool> dummyBool;
  //at least one column is viable. continue search for each viable column.
  for (auto i = 0; i != vCols.size(); ++i) {
    currentDataVector = currentDataMat[vCols[i]];
    cutPoints = getCutPointsForViableColumn(vCols[i],sortedDataMat,dThreshold);
    rCols = activeColumnsInSubset;
    if (repeatedSplittingAllowed) {
      //the user has indicated that a candidate cluster search can search a long a column vector multiple times.
      for (auto v = 0; v != rCols.size(); ++v) {
	rCols[v] = true;
      }
    }
    rCols[vCols[i]]=false; //prevent searching along same column twice (in a row for repeated splits).
    if (cutPoints.size() > maximumGateNum)
      continue;
    else {
      columnHasLowRes = lowRestrictionVector[vCols[i]];
      columnHasHighRes = highRestrictionVector[vCols[i]];
      
      //record information about intermediate gate locations for later processing.
      if (recordIndices) {
	//this is to record indices for later conditional plots in R.
	dummyBool.clear();
	for (auto b : activeRowsInSubset)
	  dummyBool.push_back(b);
      }
      else {
	dummyBool.clear();
	dummyBool.push_back(false);
      }
      currentGateInfo = {vCols[i], cutPoints.size(), cutPoints, nodeDepth, pPathScore, nodePopCounts[(vCols[i])], dummyBool};
      //std::cout << "cGateInfo: Col " << ((currentGateInfo.colNumber) + 1) << "; depth " <<
      //currentGateInfo.gateDepth << "; score " << currentGateInfo.gateDipPathScore << std::endl;;
      gdvRef.push_back(currentGateInfo);
      childPathScore = (pPathScore * (1-vColScores[i]));
      //finally, update stack with new subsets.
      for (auto j = 1; j != cutPoints.size(); ++j) {
	if ((j == 1) && (columnHasLowRes)) {
	  //the user signaled restrictrion from below. In this case, it is possible that
	  //cutPoints[0] > min(currentDataVector), since the taut string operates on the
	  //*admissible* values. In such cases, all values in currentDataVector between
	  //the minimum value and the first gate should be clusetered together.
	  gateLow = *std::min_element(currentDataVector.begin(),currentDataVector.end());
	}
	else {
	  gateLow = cutPoints[(j-1)];
	}
	if ((j == (cutPoints.size()-1)) && (columnHasHighRes)) {
	  //the user signaled restrictrion from above. In this case, it is possible that
	  //cutPoints[(cutPoints.size()-1)] < max(currentDataVector), since the taut string operates on the
	  //*admissible* values. In such cases, all values in currentDataVector between
	  //the maximum value and the last gate should be clusetered together.
	  gateHigh = *std::max_element(currentDataVector.begin(),currentDataVector.end());
	}
	else {
	  gateHigh = cutPoints[j];
	}
	rRows = activeRowsInSubset;
	//find data that falls between gates. note here we read from the unsorted active data, so that
	//we recover the observation's index in the underlying data structure.
	for (auto k = 0; k != currentDataVector.size(); ++k) {
	  if ((currentDataVector[k] <= gateLow) || (currentDataVector[k] > gateHigh)) {
	    rRows[baseIndex[k]] = false;
	  }
	  if ((j == 1) && (currentDataVector[k] == gateLow)) {
	    rRows[baseIndex[k]] = true; //the activeRow minimum of current activeCol gets assigned to the first sub-collection.
	  }
	}
	pushNode = {rRows, rCols, nodeDepth,childPathScore};
	sstackRef.push(pushNode);
      }
    }
  }
  return;
}


searchResults candidateClusterSearch(const std::vector<std::vector<double>>& dataMatrix,
				     const std::vector<std::vector<int>>& restrictionMatrix,
				     double dipThreshold,
				     int minClusterSize,
				     bool repeatedSplitting,
				     int maxDepth,
				     unsigned long long maxClusters,
				     unsigned long long maxGateNum,
				     bool randomCandSearch,
				     bool restrictedCandSearch,
				     bool annotationForestSearch,
				     unsigned long long rndSeed,
				     bool parallelEx,
				     int startingParexRoot,
				     bool recordCounts,
				     bool recordIndices) 

{
  int colNum = dataMatrix.size();
  long rowNum = (dataMatrix[0]).size();
  std::vector<bool> activeRows(rowNum,true); 
  std::vector<bool> activeCols(colNum,true); 

  std::vector<unsigned long> viableObsDenom(colNum,0);

  //initialize a small data structure to return if the counts are not recorded.
  std::vector<unsigned long> viableObsFlagVec = {0}; 
  std::vector<std::vector<unsigned long>> countViableObsMatrix(colNum,viableObsFlagVec);
  if (recordCounts) {
    //initialize data structure to store counts
    std::vector<unsigned long> viableObsFlagVecInit(rowNum,0); 
    std::vector<std::vector<unsigned long>> countViableObsMatrixInit(colNum,viableObsFlagVecInit);
    countViableObsMatrix = countViableObsMatrixInit;
    viableObsFlagVecInit.clear();
    countViableObsMatrixInit.clear();
  }

  //storage for checks on the restricted values.
  int rColCtr = 0;
  long rRowCtr; 
  
  //storage for search results
  std::vector<std::vector<long>> candidateClusters;
  std::vector<gateInfo> gateDetails;

  //intermediate search nodes are stored in the stack
  std::stack<subsetInfo> searchStack;
  struct subsetInfo currentNode;//, pushNode;
  bool firstSearchIteration = true;
  bool stillSearchingForCC = true;

  //helper variables for all search types
  int currentDepth = 0;
  double currentPathScore = 1.0;

  std::vector<bool> aRows, aCols;
  std::vector<std::vector<long>> admissibleIndices;
  std::vector<long> baseIndex, admissibleIndex;
  std::vector<std::vector<double>> activeDataSubset, admissibleDataSubset, sortedAdmissibleSubset;
  viableColumnData vcData;
  std::vector<int> viableCols;
  std::vector<double> viableColumnScores;
  std::vector<long> unresDataCounts;

  //this container used if the search is unrestricted
  unrestrictedSubsetData unresSubData;
  
  //these variables only come into use if the search has restrictions.
  std::vector<bool> lowRestrictionVec, highRestrictionVec;
  allSubsetData furSubData;

  //initialize a mersenne twister. used only if the candidate cluster search is random.
  std::mt19937 mtRCCS(rndSeed);

  //initialize flag for parallel exhuastive serch. used to forest root.
  bool firstParexIter = true;
  
  //begin candidate cluster search
  while (stillSearchingForCC) {
    if (firstSearchIteration) {
      //if it's the first iteration, all rows and columns are admissible
      firstSearchIteration = false;
      aRows = activeRows;
      aCols = activeCols;
    }
    else if (searchStack.empty()) {
      //if it's > first iteration, and the stack is empty, we're done.
      stillSearchingForCC = false;
      continue;
    }
    else {
      //grab the next node on the stack otherwise, and process it.
      currentNode = searchStack.top();
      searchStack.pop();
      aRows = currentNode.activeRowIndices;
      aCols = currentNode.activeColumnIndices;
      currentDepth = currentNode.nodeDepth;
      currentPathScore = currentNode.nodePathScore;
    }
    //if a node is too deep in a cluster tree, continue to processing next node
    if (currentDepth > maxDepth) 
      continue;
    //obtain the row indices of the observations in the active subset in the original data matrix
    baseIndex.clear();
    for (auto j = 0; j != aRows.size(); j++) 
      if (aRows[j] == true) 
	baseIndex.push_back(j);
    
    //if the current data subset is too small, continue to processing next node
    if (baseIndex.size() < minClusterSize)
      continue;
    
    if (restrictedCandSearch) {
      //if the search is restricted, in sparse data matrices it is possible for the search
      //to enter subsets where an entire entry in an active column is restricted. from
      //the perspective of the dip, this is then a null column. these
      //vacuous columns must be pruned. in the rare case that the entire subset is restricted,
      //then the search currently terminates along the branch.

      //first, examine active columns. if the number of unrestricted rows does not exceed
      //the minimum cluster size, deactivate the column
      for (auto tmpColNum = 0; tmpColNum != aCols.size(); ++tmpColNum) {
	if (aCols[tmpColNum]) {
	  rRowCtr = 0;
	  for (auto tmpRowNum = 0;  tmpRowNum != aRows.size(); ++tmpRowNum) 
	    if ((aRows[tmpRowNum]) && (((restrictionMatrix[tmpColNum])[tmpRowNum]) == 0)) 
	      rRowCtr += 1;
	  if (rRowCtr < minClusterSize) 
	    aCols[tmpColNum] = false;
	}
      }
      //now check to see if any columns remain active.
      rColCtr = 0;
      for (auto tmpColNum = 0; tmpColNum != aCols.size(); ++tmpColNum) 
	if (aCols[tmpColNum]) 
	  rColCtr += 1;
    }
    if ((restrictedCandSearch) && (rColCtr = 0)) {
      //if the 
      continue;
    }
    if (restrictedCandSearch) {
      //if the search is restricted, the *active* subset differs from the *admissible* subset.
      //nodes are updated with the *active* subset.
      //the dip and taut-string only use the *admissible* subset, corrseponding to observations
      //with restrictions==0.
      furSubData = getAllDataSubsets(aRows,aCols,dataMatrix,restrictionMatrix,recordCounts);
      activeDataSubset = furSubData.fullSubset;
      admissibleDataSubset = furSubData.urSubset;
      lowRestrictionVec = furSubData.lowResV;
      highRestrictionVec = furSubData.highResV;
      unresDataCounts = furSubData.urCountVec;
      admissibleIndices = furSubData.restrictedIndices;
    }
    else {
      //if the search is unrestricted, all observations in the active rows are admissible.
      unresSubData = getActiveDataSubset(aRows,aCols,dataMatrix);
      admissibleDataSubset = unresSubData.urSubset;
      unresDataCounts = unresSubData.urCountVec;
    }
    //sort the admissible subset for use in the dip compuations, as well as the taut-string routine.
    sortedAdmissibleSubset = sortDataSubset(admissibleDataSubset);
    //get the columns that fall below the dip threshold.
    vcData = getViableColumns(sortedAdmissibleSubset,aCols,dipThreshold,randomCandSearch,mtRCCS,
			      parallelEx,firstParexIter,startingParexRoot);
    viableCols = vcData.columnIndices;
    viableColumnScores = vcData.columnScores;
    //std::cout << "score after construction." << std::endl;
    //for (auto chk = 0; chk != viableCols.size(); ++chk) 
    // std::cout << "col " << viableCols[chk] << ": " << viableColumnScores[chk] << std::endl;
    
    //if no columns are viable, we have found a candidate cluster.
    if (viableCols.size() == 0) {
      //if not constructing an annotation forest, record it.
      if (!annotationForestSearch) {
	candidateClusters.push_back(baseIndex);
	//terminate search if we have found too many candidate cluseters
	if (candidateClusters.size() > maxClusters) {
	  stillSearchingForCC = false;
	}
      }
    }
    else {
      if (recordCounts) {
	//copy counts into ledger
	for (auto vc : viableCols) {
	  ++viableObsDenom[vc];
	  if (restrictedCandSearch) {
	    admissibleIndex = admissibleIndices[vc];
	  }
	  else {
	    admissibleIndex = baseIndex;
	  }
	  for (auto rNum : admissibleIndex) {
	    ++((countViableObsMatrix[vc])[rNum]);
	  }
	  //viableObsFlagVec = signalViableObsMatrix[vc];
	  //for (auto tv = 0; tv != viableObsFlagVec.size(); ++tv) 
	  //  if (viableObsFlagVec[tv]) 
	  //    ++((countViableObsMatrix[vc])[tv]);
	}
      }
      //update the stack, respecting the users restriction flag.
      if (restrictedCandSearch) {
	//note we pass the activeDataSubset
	restrictedStackUpdate(viableCols,aCols,aRows,activeDataSubset,sortedAdmissibleSubset,baseIndex,
			      gateDetails,searchStack,dipThreshold,maxGateNum,(currentDepth+1),
			      lowRestrictionVec,highRestrictionVec,repeatedSplitting,
			      viableColumnScores, currentPathScore, unresDataCounts,recordIndices);
      }
      else {
	//note we pass the admissibleDataSubset
	updateStackWithViableCols(viableCols,aCols,aRows,admissibleDataSubset,sortedAdmissibleSubset,baseIndex,
				  gateDetails,searchStack,dipThreshold,maxGateNum,(currentDepth+1),repeatedSplitting,
				  viableColumnScores,currentPathScore, unresDataCounts,recordIndices);
      }
    }
    if (annotationForestSearch) {
      //if we are constructing an annotation forest, termination occurs when the number of gates recorded exceeds maxClusters.
      if (gateDetails.size() > maxClusters) {
	stillSearchingForCC = false;
      }
    }
  }
  bool vacuousSearch = false;
  if (gateDetails.size() == 0) {
    // in the event no gate details are found, return the vacuous cluster (the indices of the data set)
    std::vector<long> vacuousCluster;
    for (long i = 0; i != rowNum; ++i)
      vacuousCluster.push_back(i);
    candidateClusters.push_back(vacuousCluster);
    std::vector<double> nullGates = {0.0,0.0};
    gateInfo nullInfo = {0,0,nullGates,0};
    gateDetails.push_back(nullInfo);
    if ((randomCandSearch) || (restrictedCandSearch)) {
      vacuousSearch = true;
    }
  }
  searchResults ccsResults = {candidateClusters, gateDetails, vacuousSearch, countViableObsMatrix, viableObsDenom};
  return ccsResults;
}
