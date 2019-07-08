#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <Rcpp.h>
#include <thread>
#include <mutex>
#include <condition_variable>

struct annotatedSample {
  int activeCount;
  std::vector<bool> activeRows; //flip this bit to signal a cell has been classified
  std::vector<std::string> parentAnnotation; //
  std::unordered_map<std::string,int> parentCounts;
  std::vector<std::string> childAnnotation;
  std::unordered_map<std::string,int> childCounts;
  std::unordered_map<std::string,std::unordered_set<std::string>> p2cMap;
  std::vector<std::string> activeAnnotation;
};

struct strPair {
  std::string parent;
  std::string child;
};


int getSampleCount(std::string fp) {
  //pre-processing: count how many cells are in each file
  int sampleCount = 0;
  std::string line;
  std::ifstream inputFile(fp,std::ifstream::in);
  while (getline(inputFile,line)) 
    ++sampleCount;
  inputFile.close();
  return sampleCount;
}

strPair parseAnnotations(std::string& str,
			 std::vector<int>& requestedLookupOrder,
			 std::unordered_map<std::string,std::unordered_set<std::string>>& parseMap,
			 std::vector<std::string>& requestedLookupNames)
{
  std::string parentStr="root", childStr = "";
  for (auto i = 0; i != requestedLookupOrder.size(); ++i) {
    //childStr += std::to_string(requestedLookupOrder[i]);
    childStr += requestedLookupNames[i];
    childStr += "_";
    childStr += str[requestedLookupOrder[i]];
    if (i == (requestedLookupOrder.size() - 2)) {
      parentStr = childStr;
    }
    if (i < (requestedLookupOrder.size() - 1)) {
      childStr += "_";
    }
  }
  //update the parent mapping to account for the current child
  parseMap[parentStr].insert(childStr);
  strPair outPair = {parentStr, childStr};
  return outPair;
}

bool updateAnnotationData(std::string fp,
			  std::vector<int>& requestedLookupOrder,
			  annotatedSample& sampleAnn,
			  std::vector<std::string>& requestedLookupNames)
{
  std::string line;
  std::unordered_map<std::string,std::unordered_set<std::string>> sampleMap;
  int rowNum = 0;
  int changeNum = 0;
  std::vector<std::string> sampleAnnotations(sampleAnn.activeRows.size(),"");
  std::ifstream inputFile(fp,std::ifstream::in);
  strPair tmpPair;
  while (getline(inputFile,line)) {
    if (sampleAnn.activeRows[rowNum]) {
      tmpPair = parseAnnotations(line,requestedLookupOrder,sampleMap,requestedLookupNames);
      sampleAnn.parentAnnotation[rowNum] = tmpPair.parent;
      sampleAnn.childAnnotation[rowNum] = tmpPair.child;
      ++changeNum;
    }
    ++rowNum;
  }
  inputFile.close();
  bool anyNewResults = false;
  if (changeNum) {
    anyNewResults = true;
    sampleAnn.p2cMap = sampleMap;
    int tmpCounter  = 0;
    std::string tmpStr;
    std::unordered_set<std::string> tmpChildSet;
    //update sample counts
    for (auto p : sampleAnn.p2cMap) {
      tmpCounter  = 0;
      tmpStr = p.first;
      tmpChildSet = p.second;
      for (auto i = 0 ; i != sampleAnn.parentAnnotation.size(); ++i) 
	if (sampleAnn.activeRows[i]) 
	  if (tmpStr == sampleAnn.parentAnnotation[i]) 
	    ++tmpCounter;
      sampleAnn.parentCounts[tmpStr] = tmpCounter;
      for (auto v : tmpChildSet) {
	tmpCounter = 0;
	for (auto i = 0 ; i != sampleAnn.childAnnotation.size(); ++i) 
	  if (sampleAnn.activeRows[i]) 
	    if (sampleAnn.childAnnotation[i]==v) 
	      ++tmpCounter;
	sampleAnn.childCounts[v] = tmpCounter;
      }
    }
  }
  return anyNewResults;
}

void parUpdateAnnotationData(std::string fp,
			     annotatedSample& sampleAnn,
			     int& changeCounter,
			     std::vector<int>& requestedLookupOrder,
			     std::vector<std::string>& requestedLookupNames,
			     int indexSigil,
			     std::vector<bool>& sigilVector,
			     std::condition_variable& condRef,
			     std::mutex& mutRef,
			     int& threadCounter) 
{
  //experimenting. suspect this is not thread safe
  //sampleAnn references a value in the unordered_map that represents the experiment
  //if updating sampleAnn induces a change in memory location of the annotatedSample object,
  //I think the unordered_map becomes invalidated. can eliminate this by copying the annotatedSample
  //that's being operated on to each thread, though doubles the memory footprint...
  std::string line;
  std::unordered_map<std::string,std::unordered_set<std::string>> sampleMap;
  int rowNum = 0;
  int changeNum = 0;
  std::vector<std::string> sampleAnnotations(sampleAnn.activeRows.size(),"");
  std::ifstream inputFile(fp,std::ifstream::in);
  strPair tmpPair;
  while (getline(inputFile,line)) {
    if (sampleAnn.activeRows[rowNum]) {
      tmpPair = parseAnnotations(line,requestedLookupOrder,sampleMap,requestedLookupNames);
      sampleAnn.parentAnnotation[rowNum] = tmpPair.parent;
      sampleAnn.childAnnotation[rowNum] = tmpPair.child;
      ++changeNum;
    }
    ++rowNum;
  }
  inputFile.close();
  bool anyNewResults = false;
  if (changeNum) {
    anyNewResults = true;
    sampleAnn.p2cMap = sampleMap;
    int tmpCounter  = 0;
    std::string tmpStr;
    std::unordered_set<std::string> tmpChildSet;
    //update sample counts
    for (auto p : sampleAnn.p2cMap) {
      tmpCounter  = 0;
      tmpStr = p.first;
      tmpChildSet = p.second;
      for (auto i = 0 ; i != sampleAnn.parentAnnotation.size(); ++i) 
	if (sampleAnn.activeRows[i]) 
	  if (tmpStr == sampleAnn.parentAnnotation[i]) 
	    ++tmpCounter;
      sampleAnn.parentCounts[tmpStr] = tmpCounter;
      for (auto v : tmpChildSet) {
	tmpCounter = 0;
	for (auto i = 0 ; i != sampleAnn.childAnnotation.size(); ++i) 
	  if (sampleAnn.activeRows[i]) 
	    if (sampleAnn.childAnnotation[i]==v) 
	      ++tmpCounter;
	sampleAnn.childCounts[v] = tmpCounter;
      }
    }
  }

  //once the annotatedSample is initialized, get lock for updating.
  std::lock_guard<std::mutex> guard(mutRef);
  if (anyNewResults) {
    ++changeCounter;
  }
  //update thread status, decrement the thread count, and notify scheduling thread.
  sigilVector[indexSigil] = true;
  --threadCounter;
  condRef.notify_one();
  return;
}


int tryNextPartition(const std::vector<std::string>& annotationFilePaths, 
		     const std::vector<std::string>& sampleNames,
		     std::unordered_map<std::string,annotatedSample>& sampleInfo,
		     std::vector<int> currentLookupOrder,
		     std::vector<std::string> currentLookupNames,
		     int numThreadsToUse)
{
  bool anyAnnChanges;
  int changeCount = 0;
  std::string annotationPath, sampleName;
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
    int totalFileNum = sampleNames.size();
    for (int fnum = 0; fnum != totalFileNum; ++fnum) {
      remainingFiles.push_back(fnum);
    }
    int currentFileNum;
    int usableThreadNum = std::min(totalFileNum,numThreadsToUse);
    for (int threadNum = 0; threadNum < (usableThreadNum-1); ++threadNum) {
      currentFileNum = remainingFiles.back();
      remainingFiles.pop_back();
      annotationPath = annotationFilePaths[currentFileNum];
      sampleName = sampleNames[currentFileNum];
      ++activeThreads;
      ++launchedThreads;
      threadsVec.push_back(std::thread(parUpdateAnnotationData,
				       annotationPath,
				       std::ref(sampleInfo[sampleName]),
				       std::ref(changeCount),
				       std::ref(currentLookupOrder),
				       std::ref(currentLookupNames),
				       threadNum,
				       std::ref(sigVec),
				       std::ref(parseCompleteCV),
				       std::ref(parseMutex),
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
	for (int threadNum = 0; threadNum < (usableThreadNum-1); ++threadNum) {
	  if ((sigVec[threadNum]) && (remainingFiles.size())) {
	    currentFileNum = remainingFiles.back();
	    remainingFiles.pop_back();
	    annotationPath = annotationFilePaths[currentFileNum];
	    sampleName = sampleNames[currentFileNum];
	    ++activeThreads;
	    ++launchedThreads;
	    newThreadLaunched = true;
	    sigVec[threadNum] = false;
	    if ((threadsVec[threadNum]).joinable()) {
	      (threadsVec[threadNum]).join();
	    }
	    threadsVec[threadNum] = std::thread(parUpdateAnnotationData,
						annotationPath,
						std::ref((sampleInfo[sampleName])),
						std::ref(changeCount),
						std::ref(currentLookupOrder),
						std::ref(currentLookupNames),
						threadNum,
						std::ref(sigVec),
						std::ref(parseCompleteCV),
						std::ref(parseMutex),
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
    for (auto fn = 0; fn != annotationFilePaths.size(); ++fn)  {
      sampleName = sampleNames[fn];
      annotationPath = annotationFilePaths[fn];
      anyAnnChanges = updateAnnotationData(annotationPath,currentLookupOrder,sampleInfo[sampleName],currentLookupNames);
      if (anyAnnChanges) {
	++changeCount;
      }
    }
  }
  return changeCount;
}

void getChildNodeStatus(std::unordered_map<std::string,bool>& vcMap, //update the vcMap to show if a node is vaible
			const std::unordered_set<std::string>& ccSet, //all nodes to check
			std::unordered_map<std::string,annotatedSample>& sinfoMap, //cannot be const -- unordered_map must mutate
			const std::vector<std::string>& sNames, //
			const int minSampleNum,
			const int minNodeSize){
  //clear info from previous iteration
  vcMap.clear();
  int childNodeCounter;
  annotatedSample* tmpAnn;
  std::string sampleName;
  for (auto c : ccSet) {
    childNodeCounter= 0;
    //for each child check each sample to see
    for (auto sampleNum = 0; sampleNum != sNames.size(); ++sampleNum)  {
      sampleName = sNames[sampleNum];
      tmpAnn = &(sinfoMap[sampleName]);
      if ((tmpAnn -> childCounts[c]) > minNodeSize)
	++childNodeCounter;
    }
    if (childNodeCounter >= minSampleNum) {
      vcMap[c] = true;
    }
    else {
      vcMap[c] = false;
    }
  }
  return;
}

void updateActiveAnnotation(std::unordered_map<std::string,std::unordered_set<std::string>>& gP2C,
			    std::unordered_map<std::string,bool>& vcMap, //check vcMap to see if a node is vaible
			    std::unordered_map<std::string,annotatedSample>& sinfoMap, //ref to counts
			    const std::vector<std::string>& sNames)
{
  std::unordered_set<std::string> updateParent, updateChild;
  annotatedSample* tmpAnn;
  int childCount;
  for (auto p : gP2C) { //p.first  -- parent node, p.second -- children
    childCount = 0;
    for (auto c : p.second) {
      if (vcMap[c]) {
	updateChild.insert(c);
	++childCount;
      }
    }
    //parent had no viable child nodes. 
    if (childCount == 0)
      updateParent.insert(p.first);
  }
  
  //walk the samples. update the activeAnnotation, and deactivate cells as needed.
  std::vector<std::string>* paStrVec;
  std::vector<std::string>* caStrVec;
  std::vector<std::string>* aaStrVec;
  std::vector<bool>* arVec;  
  std::string parentStr, childStr;
  int newActiveCount;
  std::string sampleName;
  for (auto sampleNum = 0; sampleNum != sNames.size(); ++sampleNum)  {
      sampleName = sNames[sampleNum];
      tmpAnn = &(sinfoMap[sampleName]);
      paStrVec = &(tmpAnn -> parentAnnotation);
      caStrVec = &(tmpAnn -> childAnnotation);
      aaStrVec = &(tmpAnn -> activeAnnotation);
      arVec = &(tmpAnn -> activeRows);
      newActiveCount = 0;
      for (auto j = 0; j != (*paStrVec).size(); ++j) {
	if (((*arVec)[j])) {
	  parentStr = (*paStrVec)[j];
	  childStr = (*caStrVec)[j];
	  if (updateParent.count(parentStr)) {
	    //have found a terminal population,record and de-activate
	    (*aaStrVec)[j] = parentStr;
	    (*arVec)[j] = false;
	  }
	  else if (updateChild.count(childStr)) {
	    //recursion continues for this pop
	    (*aaStrVec)[j] = childStr;
	    ++newActiveCount;
	  }
	  else {
	    //the cell is pruned.
	    (*aaStrVec)[j] = "0_0_0_0_0";
	    (*arVec)[j] = false;
	  }
	}
      }
      (tmpAnn -> activeCount) = newActiveCount;
  }
}


void parGetSampleCounts(int fn,
			const std::vector<std::string>& annotationFilePaths,
			const std::vector<std::string>& sampleNames,
			std::unordered_map<std::string,annotatedSample>& sampleData,
			int indexSigil,
			std::vector<bool>& sigilVector,
			std::condition_variable& condRef,
			std::mutex& mutRef,
			int& threadCounter) 
{
  //read in file, get size
  std::string annotationPath = annotationFilePaths[fn];
  int sampleSize = 0;
  std::string line;
  std::ifstream inputFile(annotationPath,std::ifstream::in);
  while (getline(inputFile,line)) 
    ++sampleSize;
  inputFile.close();
  
  //initialzie a annotatedSample structure for processing
  std::unordered_map<std::string,int> pcMap;
  std::unordered_map<std::string,int> caMap;
  std::unordered_map<std::string,std::unordered_set<std::string>> ancMap;
  std::vector<bool> arVec(sampleSize,true);
  std::vector<std::string> paVec(sampleSize,"root"); //vector of annotations per cell -- parent annotation
  std::vector<std::string> caVec(sampleSize,"root"); //vector of annotations per cell -- child annotation
  std::vector<std::string> aaVec(sampleSize,"root"); //code for a pruned cell.
  annotatedSample currentSampleInfo = {sampleSize,arVec,
				       paVec,pcMap,
				       caVec,caMap,
				       ancMap,aaVec};
  std::string sampleName = sampleNames[fn];

  //once the annotatedSample is initialized, get lock for updating.
  std::lock_guard<std::mutex> guard(mutRef);
  sampleData[sampleName] = currentSampleInfo;
  //update thread status, decrement the thread count, and notify scheduling thread.
  sigilVector[indexSigil] = true;
  --threadCounter;
  condRef.notify_one();
  return;
}


std::unordered_map<std::string,annotatedSample> initSampleContainers(const std::vector<std::string>& annotationFilePaths,
								     const std::vector<std::string>& sampleNames,
								     int numThreadsToUse)
{
  std::unordered_map<std::string,annotatedSample> sampleData; //parse container
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
    int totalFileNum = sampleNames.size();
    for (int fnum = 0; fnum != totalFileNum; ++fnum) {
      remainingFiles.push_back(fnum);
    }
    int currentFileNum;
    int usableThreadNum = std::min(totalFileNum,numThreadsToUse);
    for (int threadNum = 0; threadNum < (usableThreadNum-1); ++threadNum) {
      currentFileNum = remainingFiles.back();
      remainingFiles.pop_back();
      ++activeThreads;
      ++launchedThreads;
      threadsVec.push_back(std::thread(parGetSampleCounts,
				       currentFileNum,
				       std::ref(annotationFilePaths),
				       std::ref(sampleNames),
				       std::ref(sampleData),
				       threadNum,
				       std::ref(sigVec),
				       std::ref(parseCompleteCV),
				       std::ref(parseMutex),
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
	for (int threadNum = 0; threadNum < (usableThreadNum-1); ++threadNum) {
	  if ((sigVec[threadNum]) && (remainingFiles.size())) {
	    currentFileNum = remainingFiles.back();
	    remainingFiles.pop_back();
	    ++activeThreads;
	    ++launchedThreads;
	    newThreadLaunched = true;
	    sigVec[threadNum] = false;
	    if ((threadsVec[threadNum]).joinable()) {
	      (threadsVec[threadNum]).join();
	    }
	    threadsVec[threadNum] = std::thread(parGetSampleCounts,
						currentFileNum,
						std::ref(annotationFilePaths),
						std::ref(sampleNames),
						std::ref(sampleData),
						threadNum,
						std::ref(sigVec),
						std::ref(parseCompleteCV),
						std::ref(parseMutex),
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
    //initialize sample-level containers for single-threaded processing.
    std::unordered_map<std::string,int> pcMap;
    std::unordered_map<std::string,int> caMap;
    std::unordered_map<std::string,std::unordered_set<std::string>> ancMap;
    std::string annotationPath, sampleName;
    int sampleSize = 0;
    for (auto fn = 0; fn != annotationFilePaths.size(); ++fn)  {
      sampleName = sampleNames[fn];
      annotationPath = annotationFilePaths[fn];
      sampleSize = getSampleCount(annotationPath);
      std::vector<bool> arVec(sampleSize,true);
      std::vector<std::string> paVec(sampleSize,"root"); //vector of annotations per cell -- parent annotation
      std::vector<std::string> caVec(sampleSize,"root"); //vector of annotations per cell -- child annotation
      std::vector<std::string> aaVec(sampleSize,"root"); //code for a pruned cell.
      annotatedSample currentSampleInfo = {sampleSize,arVec,
					   paVec,pcMap,
					   caVec,caMap,
					   ancMap,aaVec};
      sampleData[sampleName] = currentSampleInfo;
    }
  }
  return sampleData;
}
								     
		   
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::IntegerMatrix getCountMatrix(std::vector<std::string> annotationFilePaths, //The first three vectors are in matching order
				   std::vector<std::string> reportingFilePaths, //sampleNames[i] is sample for reportingFilePaths[i]
				   std::vector<std::string> sampleNames, //which is path to annotationFilePaths[i]
				   std::vector<int> lookupOrder, //even valued integers corresponding to annotation cols
				   std::vector<std::string> lookupNames, //labels corresponding to annotation cols
				   int minClusterSize, //number of cells need to consider the pop for modeling
				   double minProportion,//proportion fo samples a cluster has to be in for reporting.
				   int numThreadsToUse) //parallel parsing
{
  int minSampleNumber = minProportion * annotationFilePaths.size(); //first, how many samples must a population appear in
  std::string annotationPath, sampleName;
  //create a annotatedSample structure to maintain the partition state by sample.
  //the map lookup is the sample name, represented by the string.
  std::unordered_map<std::string,annotatedSample> sampleInfo= initSampleContainers(annotationFilePaths,
										   sampleNames,
										   numThreadsToUse); 
  
  //partition the samples by the requested order
  annotatedSample* tmpAS;
  std::unordered_map<std::string,std::unordered_set<std::string>> globalP2C;
  std::unordered_map<std::string,bool> viableChildNodeMap;
  std::unordered_set<std::string> currentChildNodeSet;
  std::vector<int> currentLookupOrder;
  std::vector<std::string> currentLookupNames;
  int anyChanges = 0;
  for (auto j = 0; j != lookupOrder.size(); ++j) {
    currentLookupOrder.push_back(lookupOrder[j]);
    currentLookupNames.push_back(lookupNames[j]);
    //update annotation data for cells that are still active.
    anyChanges =  tryNextPartition(annotationFilePaths,sampleNames,sampleInfo,
				   currentLookupOrder,currentLookupNames,numThreadsToUse);
    if (anyChanges) {
      //get the child nodes of the active data
      currentChildNodeSet.clear();
      globalP2C.clear();
      for (auto sampleNumber = 0; sampleNumber != sampleNames.size(); ++sampleNumber)  {
	sampleName = sampleNames[sampleNumber];
	tmpAS = &(sampleInfo[sampleName]);
	for (auto p : (tmpAS -> p2cMap))
	  for (auto c : p.second) {
	    currentChildNodeSet.insert(c);
	    globalP2C[p.first].insert(c);
	  }
      }
      
      //parse the samples to determine which child nodes are still viable.
      getChildNodeStatus(viableChildNodeMap,currentChildNodeSet,sampleInfo,
			 sampleNames,minSampleNumber,minClusterSize);

      //update active annotations with viable children
      updateActiveAnnotation(globalP2C,viableChildNodeMap,sampleInfo,sampleNames);
    }
    else {
      break;
    }
  }
  
  //walk the samples, collect counts, and write to csv
  std::unordered_map<std::string,std::unordered_map<std::string,int>> finalCounts;
  std::unordered_map<std::string,int> finalSampleCounts, tmpMap;
  std::string currentSampleStr;
  std::ofstream outputFile;
  for (auto sampleNumber = 0; sampleNumber != sampleNames.size(); ++sampleNumber)  {
    sampleName = sampleNames[sampleNumber];
    tmpAS = &(sampleInfo[sampleName]);
    finalSampleCounts.clear();
    currentSampleStr = (reportingFilePaths[sampleNumber] + "raggedAnnotation.csv");
    outputFile.open(currentSampleStr);
    for (auto s : (tmpAS -> activeAnnotation)) {
      ++finalSampleCounts[s];
      outputFile << s << std::endl;
    }
    outputFile.close();
    finalCounts[sampleName] = finalSampleCounts;
  }
  //construct set of all populations  -- ordered for reporting.
  std::set<std::string> finalPops;
  for (auto sampleNumber = 0; sampleNumber != sampleNames.size(); ++sampleNumber)  {
    sampleName = sampleNames[sampleNumber];
    for (auto p : finalCounts[sampleName]) 
      finalPops.insert(p.first);
  }
  std::vector<std::string> rmColumnNames;
  for (auto pop : finalPops) {
    rmColumnNames.push_back(pop);
  }
  Rcpp::IntegerMatrix resultMatrix(sampleNames.size(),finalPops.size());
  Rcpp::List dimList = Rcpp::List::create(sampleNames,rmColumnNames);
  resultMatrix.attr("dimnames") = dimList;
  int currentPopNum;
  for (auto sampleNumber = 0; sampleNumber != sampleNames.size(); ++sampleNumber)  {
    sampleName = sampleNames[sampleNumber];
    tmpMap = finalCounts[sampleName];
    currentPopNum = 0;
    for (auto pop : finalPops) { 
      resultMatrix(sampleNumber,currentPopNum) = tmpMap[pop];
      ++currentPopNum;
    }
  }
  return resultMatrix;
}
