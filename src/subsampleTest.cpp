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

//structure to hold subsampling information
struct subSampledGateData {
  std::vector<double> gateLocations;
  unsigned int numGates;
};

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
std::vector<double> subsampleTest(const std::vector<double>& sortedData,
				  unsigned long subSampleSize,
				  unsigned long subSampleIterations)
{
  
  unsigned long long rndSeed = 7123541;
  std::mt19937 mtRCCS(rndSeed);
  std::mt19937& twRef = mtRCCS;
  std::vector<double> subSortedData;
  std::vector<double> subCutPoints;
  long numIter = subSampleIterations;
  unsigned int numSubCuts;
  std::vector<subSampledGateData> subSampledGates;

  //check the difference between the observed data size and the subsample size
  double dataDiff = (std::log10(sortedData.size()) - std::log10(subSampleSize));
  
  //if we are trying to subsample less than an order of magnitude difference,
  //directly shuffle the data vector
  std::vector<long> possibleIndices,indicesToShuffle,selectedIndices;
  if (dataDiff < 0.25) {
    for (auto i = 0; i != sortedData.size(); ++i)
      possibleIndices.push_back(i);
  }
  
  //on the other hand, if the difference between the subsample size and
  //observed data size exceeds an order of magnitude, randomly select
  //according to block stucture
  std::set<unsigned long> selectedIndexSet, blockIndexSet;
  std::vector<unsigned long> samplingBounds = {0};
  unsigned long stepSize = sortedData.size()/100.0;
  while (samplingBounds.back() < ((99.0/100.0)*sortedData.size()))
    samplingBounds.push_back((samplingBounds.back()+stepSize));
  samplingBounds[(samplingBounds.size()-1)] = (sortedData.size()-1);
  unsigned long blockSampleSize = subSampleSize/100.0;
  unsigned long blockStart, blockEnd;
  
  while (numIter) {
    subSortedData.clear();
    if (dataDiff < 1) {
      indicesToShuffle = possibleIndices;
      selectedIndices.clear();
      //shuffle the indices
      std::shuffle(indicesToShuffle.begin(),indicesToShuffle.end(),twRef);
      //grab the selection
      for (auto i = 0; i != subSampleSize; ++i)
	selectedIndices.push_back(indicesToShuffle[i]);
      //sort the subsample indices
      std::sort(selectedIndices.begin(),selectedIndices.end());
      //select it
      for (auto v : selectedIndices)
	subSortedData.push_back(sortedData[v]);
    }
    else {
      //randomly select the indices according to the block structure
      selectedIndexSet.clear();
      for (auto bInd = 0; bInd != (samplingBounds.size()-1); ++bInd) {
	blockIndexSet.clear();
	blockStart = samplingBounds[bInd];
	blockEnd = samplingBounds[(bInd+1)];
	std::uniform_int_distribution<> pickBlockIndex(blockStart,blockEnd);
	while (blockIndexSet.size() < blockSampleSize)
	  blockIndexSet.insert(pickBlockIndex(twRef));
	for (auto bv : blockIndexSet) 
	  selectedIndexSet.insert(bv);
      }
      for (auto v : selectedIndexSet)  {
	//Rcpp::Rcout << v << std::endl;
	subSortedData.push_back(sortedData[v]);
      }
    }
    //gate the subsampled data and record
    subCutPoints = tsGates(subSortedData,0);
    //Rcpp::Rcout << numIter << std::endl;
    numSubCuts = static_cast<unsigned int>(subCutPoints.size());
    subSampledGateData newCuts = {subCutPoints,numSubCuts};
    subSampledGates.push_back(newCuts);
    --numIter;
  }

  //figure out which subsample appears most frequently
  std::map<unsigned int, unsigned long> gateCounter;
  for (auto v : subSampledGates) 
    ++gateCounter[v.numGates];

  unsigned long maxAppearingGate = 0;
  int selGateNum = 0;
  for (auto p : gateCounter) {
    if (p.second > maxAppearingGate) {
      maxAppearingGate = p.second;
      selGateNum = p.first;
    }
  }
  
  //for the selection, get median gate by number of gates found
  double minVal = sortedData[0];
  double maxVal = sortedData[(sortedData.size()-1)];
  std::vector<double> selectedGates,cGateValues,cGateQs;
  std::vector<double> quantileVec = {0.5};
  selectedGates.push_back(minVal);
  for (auto gn = 1; gn != (selGateNum - 1); ++gn) {
    cGateValues.clear();
    cGateQs.clear();
    for (auto v : subSampledGates) 
      if (v.numGates == selGateNum)
	cGateValues.push_back(v.gateLocations[gn]);
    cGateQs = rQuantile(cGateValues,quantileVec);
    selectedGates.push_back(cGateQs[0]);
  }
  selectedGates.push_back(maxVal);
  return selectedGates;
}
