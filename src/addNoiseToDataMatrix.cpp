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

#include "faust.h"
#include <iostream>
#include <random>
#include <vector>
#include <unordered_map>
#include <thread>
#include <mutex>
#include <Rcpp.h>

double cppRnorm(double mu, double std, std::mt19937& mtRef)
{
  std::normal_distribution<> dist(mu, std);
  return dist(mtRef);
}


bool cppRbinom(double p, std::mt19937& mtRef)
{
  std::bernoulli_distribution d(p);
  return d(mtRef);
}

double cppRunif(double a, double b, std::mt19937& mtRef) {
  std::uniform_real_distribution<double> distribution(a,b);
  return distribution(mtRef);
}


std::vector<double> addDistributionNoiseToDataVector(const std::vector<double>& dataVector,
						     bool unifNoise,
						     double scaleGaussian,
						     std::mt19937& twistRef)
{
  //This function adds shape-preserving gaussian noise to dataVector.
  std::vector<double> uniqVec = dataVector;
  std::vector<double>::iterator it;
  std::sort(uniqVec.begin(),uniqVec.end());
  // must sort before creating the unique iterator
  // otherwise non-contiguous duplicates will remain in uniqVec.
  it = std::unique(uniqVec.begin(),uniqVec.end()); 
  uniqVec.resize(std::distance(uniqVec.begin(),it));
  auto UN = uniqVec.size();
  auto dUN = UN-1;
  auto DVN = dataVector.size();
  //next, count the occurance of each unique value in the data
  double obs;
  std::vector<double> resultVec(DVN,0);
  double newVal;
  if (unifNoise) {
    //add uniform noise to remove repeated values.
    if (UN == DVN) {
      //if there are no ties, return the original vector.
      resultVec = dataVector;
    }
    else {
      //first, get counts for the dataset
      std::unordered_map<double, long> countMap;
      for (auto v : dataVector) {
	++countMap[v];
      }
      //replace duplicates by uniform samples
      double lowVal, highVal;
      for (auto k = 0; k != dataVector.size(); ++k) {
	obs = dataVector[k];
	if (countMap[obs] > 1) {
	  auto lowerUniq = std::lower_bound(uniqVec.begin(), uniqVec.end(), obs);
	  auto lowerIndex = std::distance(uniqVec.begin(), lowerUniq);
	  if (lowerIndex == 0) {
	    lowVal = uniqVec[lowerIndex];
	    highVal = ((uniqVec[lowerIndex]+uniqVec[(lowerIndex+1)])/2.0);
	    newVal = cppRunif(lowVal,highVal,twistRef);
	  }
	  else if (lowerIndex == dUN) {
	    lowVal = ((uniqVec[lowerIndex]+uniqVec[(lowerIndex-1)])/2.0);
	    highVal = uniqVec[lowerIndex];
	    newVal = cppRunif(lowVal,highVal,twistRef);
	  }
	  else {
	    lowVal = ((uniqVec[lowerIndex]+uniqVec[(lowerIndex-1)])/2.0);
	    highVal = ((uniqVec[lowerIndex]+uniqVec[(lowerIndex+1)])/2.0);
	    newVal = cppRunif(lowVal,highVal,twistRef);
	  }
	}
	else {
	  newVal = obs;
	}
	resultVec[k] = newVal;
      }
    }
  }
  else {
    //adding whp order-preserving gaussian noise
    //gaussians centered at an observation.
    //sds based on distances between uniq observations.
    //record two copies to avoid later indexing errors.
    std::vector<double> dLow(UN,0.0);
    std::vector<double> dHigh(UN,0.0);
    for (auto j=0; j != UN; ++j) {
      if (j == 0) {
	dLow[j] = 0;
	dHigh[j] = uniqVec[(j+1)]-uniqVec[j];
      }
      else if (j == dUN) {
	dLow[j] = uniqVec[j]-uniqVec[(j-1)];
	dHigh[j] = 0;
      }
      else {
	dLow[j] = uniqVec[j]-uniqVec[(j-1)];
	dHigh[j] = uniqVec[(j+1)]-uniqVec[j];
      }
    }
    //sample gaussian
    double sdVal;
    for (auto k = 0; k != DVN; ++k) {
      obs = dataVector[k];
      auto lowerUniq = std::lower_bound(uniqVec.begin(), uniqVec.end(), obs);
      auto lowerIndex = std::distance(uniqVec.begin(), lowerUniq);
      if (lowerIndex == 0) {
	sdVal = (dHigh[lowerIndex]/scaleGaussian); //scaleGaussian was fixed at 4.0
      }
      else if (lowerIndex == dUN) {
	sdVal = (dLow[lowerIndex]/scaleGaussian); //scaleGaussian was fixed at 4.0
      }
      else {
	sdVal = ((dHigh[lowerIndex]+dLow[lowerIndex])/(scaleGaussian * 2.0)); //was fixed at 8.0
      }
      newVal = cppRnorm(obs,sdVal,twistRef);
      resultVec[k] = newVal;
    }
  }
  return resultVec;
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
std::vector<double> addNoiseToDataVector(const std::vector<double>& dataVector, double scaleGaussian,
					 unsigned long long rSeed) {
  std::mt19937 mt(rSeed);
  std::vector<double> dataVecWithUnifNoise = addDistributionNoiseToDataVector(dataVector, true, scaleGaussian,mt);
  std::vector<double> dataWithAllNoise = addDistributionNoiseToDataVector(dataVecWithUnifNoise, false, scaleGaussian,mt);
  return dataWithAllNoise;
}


void updateNoiseMatrix(const std::vector<std::vector<double>>& dmRef,
		       std::vector<std::vector<double>>& rcRef,
		       std::mutex& mutRef,
		       long colIndex,
		       double scaleGaussian,
		       unsigned long long rndSeed)
{
  std::vector<double> newVec = addNoiseToDataVector(dmRef[colIndex],scaleGaussian,rndSeed);
  // create a lock_guard so threads don't collide
  std::lock_guard<std::mutex> guard(mutRef); 
  for (auto j = 0; j != newVec.size(); ++j) 
    ((rcRef[colIndex])[j]) = newVec[j];
  //lock_guard destroyed on return, unlocking rcRef for other threads.
  return;
}

std::vector<std::vector<double>> addNoiseToDataMatrix(const std::vector<std::vector<double>>& dataMatrix,
						      unsigned long numThreadsInUse,
						      double scaleGaussian,
						      unsigned long long& randomSeed)
{
  auto numDataCols = dataMatrix.size();
  //unsigned long const numThreadsInUse = std::thread::hardware_concurrency();
  std::vector<std::vector<double>> resultContainer = dataMatrix;
  std::random_device rdNDM;
  unsigned long long currentSeed = randomSeed;
  if (numThreadsInUse > 2) {
    //if user requests more than two threads, parallelize the noise
    decltype(numDataCols) currentNum = 0;
    std::mutex rcMutex;
    std::vector<std::thread> threads;
    while (currentNum != numDataCols) {
      threads.clear();  //destory existing threads (if any)
      for (auto i = 0; i < (numThreadsInUse-1); ++i) {
	if (currentNum < numDataCols) {
	  // if the user has set a random seed, progressively increment for each thread.
	  if (randomSeed > 0) {
	    currentSeed += 1;
	  }
	  //otherwise, set seed for thread from random device.
	  else {
	    currentSeed = rdNDM();
	  }
	  //if there are still columns to noise, make thread to add noise to one.
	  threads.push_back(std::thread(updateNoiseMatrix,
					std::ref(dataMatrix),
					std::ref(resultContainer),
					std::ref(rcMutex),
					currentNum,
					scaleGaussian,
					currentSeed));
	  ++currentNum;
	}
	else {
	  //we have added noise to all columns
	  continue;
	}
      }
      //join all threads we spun up
      std::for_each(threads.begin(),threads.end(),
		    std::mem_fn(&std::thread::join));
    }
  }
  else {
    for (auto i = 0; i != numDataCols; ++i) {
      // if the user has set a random seed, progressively increment for each column in dataMatrix.
      if (randomSeed > 0) {
	currentSeed += 1;
      }
      //otherwise, set seed for each column from random device.
      else {
	currentSeed = rdNDM();
      }
      std::vector<double> newVec = addNoiseToDataVector(dataMatrix[i],scaleGaussian,currentSeed);
      for (auto j = 0; j != newVec.size(); ++j) {
	((resultContainer[i])[j]) = newVec[j];
      }
    }
  }
  //in either event, update the random seed if it has been incremented.
  if (randomSeed > 0) {
    randomSeed = currentSeed;
  }
  return resultContainer;
}
