#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <Rcpp.h>
#include <thread>
#include <mutex>
#include <condition_variable>

std::string parseAnnotations(std::string& str,
			     std::vector<int>& requestedLookup) {
  std::string strOut = "";
  for (auto i = 0; i != requestedLookup.size(); ++i) {
    strOut += str[requestedLookup[i]];
    if (i < (requestedLookup.size() - 1)) {
      strOut += "_";
    }
  }
  return strOut;
}

std::unordered_map<std::string,int> getSampleAnns(std::string fp,
						  std::vector<int>& requestedLookup) {
  std::string line;
  std::unordered_map<std::string,int> sampleAnns;
  std::ifstream inputFile(fp,std::ifstream::in);
  while (getline(inputFile,line)) 
    ++sampleAnns[parseAnnotations(line,requestedLookup)];
  inputFile.close();
  return sampleAnns;
}


void  parGetSampleAnns(std::string fp,
		       std::vector<int> requestedLookup,
		       int indexSigil,
		       std::vector<bool>& sigilVector,
		       std::condition_variable& condRef,
		       std::mutex& mutRef,
		       std::vector<std::unordered_map<std::string,int>>& parseContainer,
		       int& threadCounter) 
{
  std::string line;
  std::unordered_map<std::string,int> sampleAnns;
  std::ifstream inputFile(fp,std::ifstream::in);
  while (getline(inputFile,line)) { 
    ++sampleAnns[parseAnnotations(line,requestedLookup)];
  }
  inputFile.close();
  //parse complete. get lock for updating.
  std::lock_guard<std::mutex> guard(mutRef); 
  parseContainer.push_back(sampleAnns);
  sigilVector[indexSigil] = true;
  //decrement the thread count, and notify scheduling thread.
  --threadCounter;
  condRef.notify_one();
  return;
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::List getAnnotationStats(std::vector<std::string> filePaths, //paths to the annotation matrix csv
			      std::vector<int> requestedLookup, //even valued integers corresponding to annotation cols
			      int minClusterSize, //number of cells need to consider the pop for modeling
			      double minProportion,//proportion fo samples a cluster has to be in for reporting.
			      int numThreadsToUse) //parallel parsing
{
  int minSampleNumber = minProportion * filePaths.size();
  std::vector<std::unordered_map<std::string,int>> expAnns;
  std::unordered_map<std::string,int> sampAnns, expSummary;
  //parse each sample's annotation matrix
  if (numThreadsToUse > 1) {
    //parallel parse. bring in locks and notification machinery to control threads.
    std::mutex parseMutex; 
    std::unique_lock<std::mutex> guardParseUpdate(parseMutex); 
    std::condition_variable parseCompleteCV;
    int launchedThreads = 0, activeThreads = 0;
    std::vector<bool> sigVec((numThreadsToUse-1),false);
    std::vector<std::thread> threadsVec;
    bool newThreadLaunched = false;
    std::vector<int> remainingFiles;
    int totalFileNum = filePaths.size();
    for (int fnum = 0; fnum != totalFileNum; ++fnum) {
      remainingFiles.push_back(fnum);
    }
    int currentFileNum;
    std::string currentFilePath;
    int usableThreadNum = std::min(totalFileNum,numThreadsToUse);
    for (int fileNum = 0; fileNum < (usableThreadNum-1); ++fileNum) {
      currentFileNum = remainingFiles.back();
      currentFilePath = filePaths[currentFileNum];
      remainingFiles.pop_back();
      ++activeThreads;
      ++launchedThreads;
      threadsVec.push_back(std::thread(parGetSampleAnns,
				       currentFilePath,
				       requestedLookup,
				       fileNum,
				       std::ref(sigVec),
				       std::ref(parseCompleteCV),
				       std::ref(parseMutex),
				       std::ref(expAnns),
				       std::ref(activeThreads)));
					 
    }
    if (totalFileNum < numThreadsToUse) {
      while (launchedThreads) {
	parseCompleteCV.wait(guardParseUpdate,
			     [&](){return !(activeThreads == launchedThreads);});
	launchedThreads = activeThreads; //each thread terminates by decrementing activeTrhreads.
      }
    }
    else {
      //monitor thread pool, launching new threads along new roots as each search completes.
      while (launchedThreads) {
	parseCompleteCV.wait(guardParseUpdate,
			     [&](){return !(activeThreads == launchedThreads);});
	launchedThreads = activeThreads;
	newThreadLaunched = false;
	for (int fileNum = 0; fileNum < (usableThreadNum-1); ++fileNum) {
	  if ((sigVec[fileNum]) && (remainingFiles.size())) {
	    currentFileNum = remainingFiles.back();
	    currentFilePath = filePaths[currentFileNum];
	    remainingFiles.pop_back();
	    ++activeThreads;
	    ++launchedThreads;
	    newThreadLaunched = true;
	    sigVec[fileNum] = false;
	    if ((threadsVec[fileNum]).joinable()) {
	      (threadsVec[fileNum]).join();
	    }
	    threadsVec[fileNum] = std::thread(parGetSampleAnns,
					      currentFilePath,
					      requestedLookup,
					      fileNum,
					      std::ref(sigVec),
					      std::ref(parseCompleteCV),
					      std::ref(parseMutex),
					      std::ref(expAnns),
					      std::ref(activeThreads));
	  }
	}
      }
      for (int i = 0; i < (usableThreadNum-1); ++i) {
	if ((threadsVec[i]).joinable()) {
	  (threadsVec[i]).join();
	}
      }
    }
  }
  else {
    for (auto fn = 0; fn != filePaths.size(); ++fn)  {
      expAnns.push_back(getSampleAnns(filePaths[fn],requestedLookup));
    }
  }
  
  //examine the experimental summary for large subsets
  for (auto sn = 0; sn != expAnns.size(); ++sn) {
    for (auto p : expAnns[sn]) {
      if (p.second > minClusterSize) {
	++expSummary[p.first];
      }
    }
  }
  
  //accumulate clusters that appear in more than minSampleNumber samples.
  std::vector<int> clusterAppNumVec;
  std::vector<std::string> foundClusters;
  for (auto p : expSummary) {
    if (p.second >= minSampleNumber) {
      clusterAppNumVec.push_back(p.second);
      foundClusters.push_back(p.first);
    }
  }
  
  if (foundClusters.size() == 0) {
    clusterAppNumVec.push_back(0);
    foundClusters.push_back("emptyPartition");
  }
  Rcpp::List outList = Rcpp::List::create(Rcpp::Named("Queries")=foundClusters,
					  Rcpp::Named("numSamples")=clusterAppNumVec);
  return outList;
    //return Rcpp::wrap(foundClusters);
}
