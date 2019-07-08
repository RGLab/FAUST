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
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <Rcpp.h>
#include "faust.h"
#include <sys/resource.h>
#include <thread>

#ifndef _DEBUG_USER_ANNOTATION_FOREST_
#define _DEBUG_USER_ANNOTATION_FOREST_ false
#endif

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::List cppGrowAnnotationForest(Rcpp::NumericMatrix& rawDataMatrix,
				   double dipT,
				   int clusterLB,
				   bool repeatsAllowed,
				   int maxSearchDepth,
				   unsigned long long maxClusterNum,
				   bool verboseForestRun,
				   unsigned long long maxNumberOfGates,
				   bool randomSearch,
				   int numThreadsRequested,
				   bool useRestrictedValue,
				   Rcpp::NumericMatrix& restrictedValueMatrix,
				   int numGateUB,
				   double maxSearchTime,
				   double gaussianScale,
				   unsigned long long randomSeed,
				   bool recordCounts,
				   bool recordIndices)

{
  //Determine number of threads to use.
  unsigned long numThreadsToUse;
  unsigned long const maxNumThreadsPossible = std::thread::hardware_concurrency();
  if (numThreadsRequested < 1) {
    if (verboseForestRun) {
      Rcpp::Rcout << "User requested max threads." << std::endl;
      Rcpp::Rcout << "Hardware supports maximum of " << maxNumThreadsPossible << " threads." << std::endl;
    }
    numThreadsToUse = maxNumThreadsPossible;
    if (verboseForestRun) {
      Rcpp::Rcout << "Using " << numThreadsToUse << " threads." << std::endl;
    }
  }
  else {
    if (verboseForestRun) {
      Rcpp::Rcout << "User requested " << numThreadsRequested << ". " << std::endl;
      Rcpp::Rcout << "Hardware supports maximum of " << maxNumThreadsPossible << " threads." << std::endl;
    }
    numThreadsToUse = numThreadsRequested;
    if (verboseForestRun) {
      Rcpp::Rcout << "Using " << numThreadsToUse << " threads." << std::endl;
    }
  }

  if (verboseForestRun) {
    Rcpp::Rcout << "Copying data into C++ data structures." << std::endl;
  }
  //To begin analysis, copy the R Matrix into desired C++ data structure
  std::vector<std::vector<double>> dataValsIn;
  std::vector<std::vector<int>> restrictedVals;
      
  std::vector<double> tmpV;
  std::vector<int> tmpV2;
  auto colNum = rawDataMatrix.ncol();
  for (auto i = 0; i != colNum; ++i) {
    tmpV.clear();
    tmpV2.clear();
    Rcpp::NumericMatrix::Column tmp = rawDataMatrix(Rcpp::_,i);
    for (auto j : tmp)
      tmpV.push_back(j);
    dataValsIn.push_back(tmpV);
    if (useRestrictedValue) {
      Rcpp::NumericMatrix::Column tmp2 = restrictedValueMatrix(Rcpp::_,i);
      for (auto k : tmp2) {
	if (k == 1) {
	  tmpV2.push_back(1);
	}
	else if (k ==2) {
	  tmpV2.push_back(2);
	}
	else {
	  //default to unrestricted in the event a non-{0,1,2} int is in set.
	  //R wrapper is responsible for making sure this does not happen.
	  tmpV2.push_back(0);
	}
      }
      restrictedVals.push_back(tmpV2);
    }
  }
  if (verboseForestRun) {
    Rcpp::Rcout << "Copy complete." << std::endl;
    Rcpp::Rcout << "Beginning to add noise." << std::endl;
  }
  unsigned long long currentRandomSeed = randomSeed;
  //Next, add Uniform and Gaussian noise to input data.
  std::vector<std::vector<double>> dataVals = addNoiseToDataMatrix(dataValsIn,numThreadsToUse,gaussianScale,
								   currentRandomSeed);
  if (verboseForestRun) {
    Rcpp::Rcout << "Noise added." << std::endl;
    Rcpp::Rcout << "Growing annotation forest." << std::endl;
  }
  bool runningAnnotationForest = true;
  //Next, we search for candidate clusters.
  searchResults candClusters = findCandidateClusters(dataVals,
						     restrictedVals,
						     dipT,
						     clusterLB,
						     repeatsAllowed,
						     maxSearchDepth,
						     maxClusterNum,
						     maxNumberOfGates,
						     useRestrictedValue,
						     randomSearch,
						     maxSearchTime,
						     verboseForestRun,
						     runningAnnotationForest,
						     numThreadsToUse,
						     currentRandomSeed,
						     recordCounts,
						     recordIndices);

  if (verboseForestRun) {
    Rcpp::Rcout << "Annotation forest grown." << std::endl;
    Rcpp::Rcout << "Parsing and returning annotation forest." << std::endl;
  }
  
  //record gate locations found by the candidate cluster search
  std::vector<gateInfo> gatePlacements = candClusters.gateLocations;
  std::vector<std::vector<double>> gatesByColumn;
  gatesByColumn.resize((5*colNum));
  std::vector<std::vector<bool>> indicesByColumn;
  indicesByColumn.resize(colNum);  
  //double count the depths so when the indices are recorded, we can pull out the correct depths.
  std::vector<std::vector<int>> depthsByColumn;
  depthsByColumn.resize(colNum);  
  gateInfo currentGateInfo;
  int gateColumnNum, indexColumnNum;
  int currentGateDepth;
  unsigned long currentNumGates;
  double currentNodeScore;
  long currentNodePopSize;
  std::vector<double> gateLocs;
  std::vector<bool> gateIndices;
  for (auto i = 0; i != gatePlacements.size(); ++i) {
    currentGateInfo = gatePlacements[i];
    gateColumnNum = (5*(currentGateInfo.colNumber));
    gateLocs = currentGateInfo.gates;
    gateIndices = currentGateInfo.gateActiveRows;
    indexColumnNum = currentGateInfo.colNumber;
    currentNumGates = ((currentGateInfo.numGates)-2);
    currentGateDepth = currentGateInfo.gateDepth;
    //always record indices for conditional plotting.
    //if recordIndices = false, then the cost is one bit per currentGateInfo
    //we do incur the cost of returning a list of ints for depthsByColumn
    depthsByColumn[indexColumnNum].push_back(currentGateDepth);
    for (auto b : gateIndices) 
      indicesByColumn[indexColumnNum].push_back(b);
    //Rcpp::Rcout << "node depth: " << currentGateDepth << "; node score: " << currentGateInfo.gateDipPathScore <<
    //  "; node size: " << currentGateInfo.nodeObsCount << std::endl;
    currentNodeScore = currentGateInfo.gateDipPathScore;
    currentNodePopSize = currentGateInfo.nodeObsCount;
    if (currentNumGates > numGateUB) {
      continue; //user supplies the upper bound to ignore.
    }
    else {
      for (auto j = 1; j != (gateLocs.size()-1); ++j) {
	(gatesByColumn[gateColumnNum]).push_back(gateLocs[j]);
	(gatesByColumn[(gateColumnNum+1)]).push_back(currentNumGates);
	(gatesByColumn[(gateColumnNum+2)]).push_back(currentGateDepth);
	(gatesByColumn[(gateColumnNum+3)]).push_back(currentNodeScore);
	(gatesByColumn[(gateColumnNum+4)]).push_back(currentNodePopSize);
      }
    }
  }
  Rcpp::List finalGateData = Rcpp::wrap(gatesByColumn);
  Rcpp::List finalIndexData = Rcpp::wrap(indicesByColumn);
  Rcpp::List finalIndexDepthData = Rcpp::wrap(depthsByColumn);
  std::vector<std::vector<unsigned long>> finalCounts = candClusters.subsetCounts;
  std::vector<unsigned long> tmpFC;

  //Rcpp::Rcout << "final counts stats" <<std::endl;
  //Rcpp::Rcout << "cols: "<< finalCounts.size() << std::endl;
  //for (auto chk = 0; chk != finalCounts.size(); ++chk)
  // Rcpp::Rcout << "row " << chk << ": " << finalCounts[chk].size() << std::endl;

  //default matrix to return
  Rcpp::NumericMatrix finalCountMatrix(2,colNum);
  if (recordCounts) {
    //get count data from annotation forest
    Rcpp::NumericMatrix finalCountMatrixRecord(((finalCounts[0]).size()),colNum);
    for (auto tcNum = 0; tcNum != colNum; ++tcNum) {
      tmpFC = finalCounts[tcNum];
      for (auto tmpRN = 0; tmpRN != tmpFC.size(); ++tmpRN)
	finalCountMatrixRecord(tmpRN,tcNum) = tmpFC[tmpRN];
    }
    finalCountMatrix = finalCountMatrixRecord;
  }

  Rcpp::NumericVector finalSubsetDenom = Rcpp::wrap(candClusters.subsetDenom);
  Rcpp::List outList = Rcpp::List::create(Rcpp::Named("gateData")=finalGateData,
					  Rcpp::Named("subsetCounts")=finalCountMatrix,
					  Rcpp::Named("subsetDenom")=finalSubsetDenom,
					  Rcpp::Named("indexData")=finalIndexData,
					  Rcpp::Named("indexDepthData")=finalIndexDepthData);
  return outList;
}
