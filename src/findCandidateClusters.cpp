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
#include <cmath>
#include <algorithm>
#include "faust.h"
#include <chrono>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <sys/resource.h>
#include <random>

void ccSearchThread(const std::vector<std::vector<double>>& dmRef,
		    const std::vector<std::vector<int>>& rvRef,
		    const double& dipVal,
		    const int& cLowerBound,
		    const bool& repAllow,
		    const int& maxSearchD,
		    const unsigned long long& mclNum,
		    const unsigned long long& maxNumGat,
		    std::mutex& mutRef,
		    searchResults& ccdsRef,
		    const bool& searchRandom,
		    const bool& searchRestricted,
		    long& bsCounter,
		    const bool& annotationForest,
		    const long& pathologyLimit,
		    std::vector<bool>& sigilVector,
		    const int indexSigil,//important to copy this index from main thread -- else iteration points past end of sigilVector
		    unsigned long long rSeed,
		    const bool& searchParEx,
		    int parExStartRoot,//also need to copy in, otherwise roots repeat.
		    std::condition_variable& condRef,
		    int& threadCounter,
		    bool recordCounts,
		    bool recordIndices)

{
  //randomly sample from the space of clustering trees, subject to the restrictions in rvRef
  searchResults parResult =  candidateClusterSearch(dmRef,
						    rvRef,
						    dipVal,
						    cLowerBound,
						    repAllow,
						    maxSearchD,
						    mclNum,
						    maxNumGat,
						    searchRandom,
						    searchRestricted,
						    annotationForest,
						    rSeed,
						    searchParEx,
						    parExStartRoot,
						    recordCounts,
						    recordIndices);
  
  //having completed a random candidate cluster search, lock for updating
  std::lock_guard<std::mutex> guard(mutRef); 
  if ((parResult.abortIteration) && (((dmRef[0]).size()) > (100*pathologyLimit))) {
    //this condition is only true if the search didn't find a candidate cluster
    //and the subset size exceed a multiple of the pathology limit. The size check is b/c
    //we want to allow residual clusters to be lumped together.
    ++bsCounter;
  }
  else {
    std::vector<std::vector<unsigned long>> subCtr = parResult.subsetCounts;
    std::vector<unsigned long> tmpCol;
    if (recordCounts) {
      for (auto tcNum = 0; tcNum != subCtr.size(); ++tcNum) {
	tmpCol = subCtr[tcNum];
	for (auto trNum = 0; trNum != tmpCol.size(); ++trNum)
	  (((ccdsRef.subsetCounts)[tcNum])[trNum]) += tmpCol[trNum];
      }
      tmpCol = parResult.subsetDenom;
      for (auto tcNum = 0; tcNum != tmpCol.size(); ++tcNum) {
	((ccdsRef.subsetDenom)[tcNum]) += tmpCol[tcNum];
      }
    }
    //Rcpp::Rcout << "update cc" << std::endl;
    std::vector<std::vector<long>> parCands = parResult.candidates;
    for (auto cc : parCands) {
      ccdsRef.candidates.push_back(cc);
    }
    //Rcpp::Rcout << "update gl" << std::endl;
    std::vector<gateInfo> parGateLocs = parResult.gateLocations;
    //Rcpp::Rcout << "gate info found in thread " << indexSigil << std::endl;
    for (auto gL : parGateLocs) {
      ccdsRef.gateLocations.push_back(gL);
      //Rcpp::Rcout << "colNumber: " << gL.colNumber << "; numGates: " << gL.numGates << "; gateDepth: " << gL.gateDepth << std::endl;
    }
    ccdsRef.abortIteration = false;
  }
  //update the sigil vector to notify main thread this thread can be relaunched.
  sigilVector[indexSigil] = true;
  //decrement the thread count, and notify scheduling thread.
  --threadCounter;
  condRef.notify_one();
  return;
}


searchResults findCandidateClusters(const std::vector<std::vector<double>>& fDataVals,
				    const std::vector<std::vector<int>>& fRestrictedVals,
				    const double& fDipT,
				    const int& fClusterLB,
				    const bool& fRepeatsAllowed,
				    const int& fMaxSearchDepth,
				    const unsigned long long& fMaxClusterNum,
				    const unsigned long long& fMaxNumberOfGates,
				    const bool& fUseRestrictedValue,
				    const bool& fRandomSearch,
				    const double& maxAllowedTime,
				    const bool& printDebugInfo,
				    const bool& annotationForestRun,
				    unsigned long numThreadsTotal,
				    unsigned long long& randomSeed,
				    bool recordCounts,
				    bool recordIndices)
{
  searchResults fCandClusters;
  int fColNum = fDataVals.size();
  long fRowNum = (fDataVals[0]).size();
  std::vector<unsigned long> initSubsetDenom(fColNum,0);

  //small data structure initialized -- recordCounts == false
  std::vector<unsigned long> initSubsetCounter = {0};
  std::vector<std::vector<unsigned long>> searchSubsetCounter(fColNum,initSubsetCounter); 	
  if (recordCounts) {
    std::vector<unsigned long> initSubsetCounterAll(fRowNum,0);
    std::vector<std::vector<unsigned long>> searchSubsetCounterInit(fColNum,initSubsetCounterAll);
    searchSubsetCounter = searchSubsetCounterInit;
    searchSubsetCounterInit.clear();
    initSubsetCounterAll.clear();
  }
  fCandClusters.subsetCounts = searchSubsetCounter;
  fCandClusters.subsetDenom = initSubsetDenom;
  bool underMaxTime = true;
  bool noPathology = true;
  auto startSearch = std::chrono::steady_clock::now();
  auto checkPoint = std::chrono::steady_clock::now();
  std::chrono::duration<double> currentDuration;
  double elapsedTime;
  long badSearchCounter = 0;
  //const auto numInSubset = ((fDataVals[0]).size());
  const long pathologyLimit = (fClusterLB * 10); //only check for vacuous clusterings if the subset is large enough.
  searchResults parResult;
  bool noAnnForestInterrupt = true;
  //auto duration_s = std::chrono::duration<double>(dcast);
  if (printDebugInfo) {
    Rcpp::Rcout << "Max allowed time for search: " << maxAllowedTime << std::endl;
  }
  std::random_device rdCCS;
  unsigned long long currentSeed = randomSeed;
  bool ccSearchParex = false; // flag for exhaustive parallel search.
  int ccSearchParexRoot = 0; //used to signal starting root for exhaustive parallel search.
  if (numThreadsTotal > 2) {
    std::mutex ccMutex;
    std::unique_lock<std::mutex> guardCandidateClusters(ccMutex);
    std::condition_variable searchCompleteCV;
    int launchedThreads = 0, activeThreads = 0;
    std::vector<bool> sigVec((numThreadsTotal-1),false);
    std::vector<std::thread> threadsVec;
    bool newThreadLaunched = false;
    //build up vector of starting roots.
    std::vector<bool> vRoots = determineAnnotationColumns(fDataVals,fDataVals.size(),
							  fDipT,fRestrictedVals,
							  fUseRestrictedValue);
    std::vector<int> vRootsIndex;
    for (int vri=0; vri !=vRoots.size();++vri)
      if (vRoots[vri])
	vRootsIndex.push_back(vri);
    if (printDebugInfo) {
      Rcpp::Rcout << "Searching across " << vRootsIndex.size() << " roots in forest." << std::endl;
    }
    int allRootNum = vRootsIndex.size();
    if (allRootNum <= 1) {
      guardCandidateClusters.unlock();
      //in the case of no viable or a single viable root, pass to exhaustiveSearch.
      fCandClusters =  candidateClusterSearch(fDataVals,fRestrictedVals,fDipT,fClusterLB,
					      fRepeatsAllowed,fMaxSearchDepth,
					      fMaxClusterNum,fMaxNumberOfGates,fRandomSearch,fUseRestrictedValue,
					      annotationForestRun,(randomSeed+1),ccSearchParex,ccSearchParexRoot,
					      recordCounts,recordIndices);
    }
    else {
      ccSearchParex = true;
      //check if there are fewer roots in the forest than threads available
      int usableThreadNum = std::min((allRootNum+1),static_cast<int>(numThreadsTotal));
      for (int i = 0; i < (usableThreadNum-1); ++i) {
	// if the user has set a random seed, progressively increment for each thread.
	if (randomSeed > 0) {
	  currentSeed += 1;
	}
	//otherwise, set seed for thread from random device.
	else {
	  currentSeed = rdCCS();
	}
	ccSearchParexRoot = vRootsIndex.back();
	vRootsIndex.pop_back();
	++activeThreads;
	++launchedThreads;
	if (printDebugInfo) {
	  Rcpp::Rcout << "fRandomSearch at launch: " << fRandomSearch << std::endl;
	  Rcpp::Rcout << "ccSearchParex at launch: " << ccSearchParex << std::endl;
	  Rcpp::Rcout << "ccSearchParexRoot at launch: " << ccSearchParexRoot << "." << std::endl;
	  Rcpp::Rcout << "launchedThreads at launch: " << launchedThreads << "." << std::endl;
	  Rcpp::Rcout << "activeThreads at launch: " << activeThreads << "." << std::endl;
	}
	threadsVec.push_back(std::thread(ccSearchThread,
					 std::ref(fDataVals),
					 std::ref(fRestrictedVals),
					 std::ref(fDipT),
					 std::ref(fClusterLB),
					 std::ref(fRepeatsAllowed),
					 std::ref(fMaxSearchDepth),
					 std::ref(fMaxClusterNum),
					 std::ref(fMaxNumberOfGates),
					 std::ref(ccMutex),
					 std::ref(fCandClusters),
					 std::ref(fRandomSearch),
					 std::ref(fUseRestrictedValue),
					 std::ref(badSearchCounter),
					 std::ref(annotationForestRun),
					 std::ref(pathologyLimit),
					 std::ref(sigVec),
					 i,
					 currentSeed,//(randomSeed+1),
					 std::ref(ccSearchParex),
					 ccSearchParexRoot,
					 std::ref(searchCompleteCV),
					 std::ref(activeThreads),
					 recordCounts,
					 recordIndices));
	if (printDebugInfo) {
	  Rcpp::Rcout << "launched thread " << i << std::endl;
	  Rcpp::Rcout << "Used seed " << currentSeed << std::endl;//(randomSeed+1) << std::endl;
	}
      }
      if (vRootsIndex.size() == 0) {
	if (printDebugInfo) {
	  Rcpp::Rcout << "More theads than roots. Spin off." << std::endl;
	  Rcpp::Rcout << "launchedThreads at start: " << launchedThreads << "." << std::endl;
	  Rcpp::Rcout << "activeThreads at start: " << activeThreads << "." << std::endl;
	}
	//we have launched searches along the entire space. wait for them to complete.
	while (launchedThreads) {
	  searchCompleteCV.wait(guardCandidateClusters,
				[&](){return !(activeThreads == launchedThreads);});
	  if (printDebugInfo) {
	    Rcpp::Rcout << " after alert, launchedThreads: " << launchedThreads << std::endl;
	    Rcpp::Rcout << " after alert, activeThreads: " << activeThreads << std::endl;
	  }
	  launchedThreads = activeThreads;
	  if (printDebugInfo) {
	    Rcpp::Rcout << " after reconcile, launchedThreads: " << launchedThreads << std::endl;
	    Rcpp::Rcout << " after reconcile, activeThreads: " << activeThreads << std::endl;
	  }
	}
	if (printDebugInfo)
	  Rcpp::Rcout << "search complete: " << launchedThreads << " remain." << std::endl;
      }
      else {
	//monitor thread pool, launching new threads along new roots as each search completes.
	while (launchedThreads && underMaxTime && noPathology && noAnnForestInterrupt) {
	  //as in the previous case, each thread terminates by decrementing activeThreads.
	  //here, by starting a loop iteration with a wait, and then setting launchedThreads to activeThreads
	  //loop will exit once activeThreads have explored the remainder of vRootsIndex.
	  searchCompleteCV.wait(guardCandidateClusters,
				[&](){return !(activeThreads == launchedThreads);});
	  launchedThreads = activeThreads;
	  newThreadLaunched = false;
	  for (int i = 0; i < (usableThreadNum-1); ++i) {
	    if ((sigVec[i]) && (vRootsIndex.size())) {
	      //thread i has signaled completion, and there are still unexplored
	      //roots in the forest. launch another thread
	      ccSearchParexRoot = vRootsIndex.back();
	      vRootsIndex.pop_back();
	      // if the user has set a random seed, progressively increment for next thread.
	      if (randomSeed > 0) {
		currentSeed += 1;
	      }
	      //otherwise, set seed for next thread from random device.
	      else {
		currentSeed = rdCCS();
	      }
	      ++activeThreads;
	      ++launchedThreads;
	      //reset sigil and launch new thread.
	      newThreadLaunched = true;
	      sigVec[i] = false;
	      //check if the thread is joinable. if so, join it.
	      if ((threadsVec[i]).joinable()) {
		(threadsVec[i]).join();
	      }
	      //launch new thread.
	      threadsVec[i] = std::thread(ccSearchThread,
					  std::ref(fDataVals),
					  std::ref(fRestrictedVals),
					  std::ref(fDipT),
					  std::ref(fClusterLB),
					  std::ref(fRepeatsAllowed),
					  std::ref(fMaxSearchDepth),
					  std::ref(fMaxClusterNum),
					  std::ref(fMaxNumberOfGates),
					  std::ref(ccMutex),
					  std::ref(fCandClusters),
					  std::ref(fRandomSearch),
					  std::ref(fUseRestrictedValue),
					  std::ref(badSearchCounter),
					  std::ref(annotationForestRun),
					  std::ref(pathologyLimit),
					  std::ref(sigVec),
					  i,
					  currentSeed,//(randomSeed+1),
					  std::ref(ccSearchParex),
					  ccSearchParexRoot,
					  std::ref(searchCompleteCV),
					  std::ref(activeThreads),
					  recordCounts,
					  recordIndices);
	      if (printDebugInfo) {
		Rcpp::Rcout << "Launched a new thread in slot: " << i << std::endl;
		Rcpp::Rcout << "Used seed: " << currentSeed << std::endl;//(randomSeed+1) << std::endl;
	      }
	    }
	  }
	  //check if we've had too many bad searches.
	  //this check targets residual searches: numInSubset small.
	  if (badSearchCounter > (9*(fMaxClusterNum/10))) {
	    if (printDebugInfo) {
	      Rcpp::Rcout << "Search pathology: vacuous results. Abort and re-noise." << std::endl;
	    }
	    noPathology = false;
	  }
	  if (newThreadLaunched) {
	    checkPoint = std::chrono::steady_clock::now();
	    currentDuration = checkPoint-startSearch;
	    elapsedTime = currentDuration.count();
	    if (printDebugInfo) {
	      Rcpp::Rcout << "Elapsed search time: " << elapsedTime << std::endl;
	      if (annotationForestRun) {
		Rcpp::Rcout << "Gates found: " << fCandClusters.gateLocations.size() << std::endl;
	      }
	      else{
		Rcpp::Rcout << "Candidate clusters found: " << fCandClusters.candidates.size() << std::endl;
	      }
	    }
	    if (elapsedTime > maxAllowedTime) {
	      underMaxTime = false;
	    }
	    if (annotationForestRun) {
	      if (fCandClusters.gateLocations.size() > fMaxClusterNum) {
		noAnnForestInterrupt = false;
	      }
	    }
	  }
	}
      }
      //join any dangling threads
      for (int i = 0; i < (usableThreadNum-1); ++i) {
	if ((threadsVec[i]).joinable()) {
	  (threadsVec[i]).join();
	}
      }
      //get final stats
      if (printDebugInfo) {
	checkPoint = std::chrono::steady_clock::now();
	currentDuration = checkPoint-startSearch;
	elapsedTime = currentDuration.count();
	Rcpp::Rcout << "Final search time: " << elapsedTime << std::endl;
	if (annotationForestRun) {
	  Rcpp::Rcout << "Total gates found: " << fCandClusters.gateLocations.size() << std::endl;
	}
	else{
	  Rcpp::Rcout << "Total candidate clusters found: " << fCandClusters.candidates.size() << std::endl;
	}
      }
      //finally, void existing results if a search pathology or time violation is signaled.
      if ((!noPathology) || (!underMaxTime)) {
	fCandClusters.candidates.clear();
	fCandClusters.gateLocations.clear();
	fCandClusters.abortIteration = true;
      }
    }
  }
  else {
    //pass in the random seed directly in this event -- it won't be used.
    fCandClusters =  candidateClusterSearch(fDataVals,fRestrictedVals,fDipT,fClusterLB,
					    fRepeatsAllowed,fMaxSearchDepth,
					    fMaxClusterNum,fMaxNumberOfGates,fRandomSearch,fUseRestrictedValue,
					    annotationForestRun,(randomSeed+1),ccSearchParex,ccSearchParexRoot,
					    recordCounts,
					    recordIndices);

  }
  //whatever the search parameter, update the random seed with the current seed.
  if (randomSeed > 0) {
    randomSeed = currentSeed;
  }
  return fCandClusters;
}
