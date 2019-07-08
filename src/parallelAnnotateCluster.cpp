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
#include <thread>
#include <mutex>

void threadAnnotateCluster(const std::vector<std::vector<double>>& dataValsIn,
			   const std::vector<int>& clusteringIn,
			   const std::vector<bool>& annotatingCols,
			   const std::vector<std::vector<double>>& annotationBoundaries,
			   const std::vector<std::string>& annotationsIn,
			   const std::vector<double>& annComparisonQuantiles,
			   const std::vector<int>& uniqClusterLabels,
			   unsigned long clusterNum,
			   std::mutex& annMutRef,
			   std::vector<std::vector<int>>& globalAnnotatedClusters)
{
  auto numberOfColumns = dataValsIn.size();
  auto numberOfObservations = clusteringIn.size();
  std::vector<int> labelTemplatePar(numberOfColumns,0);
  std::vector<std::vector<int>> annotatedClusters(numberOfObservations,labelTemplatePar);
  //inCluster flags if an observation of dataValsIn belongs to the cluster. reset at each loop iteration.
  std::vector<bool> inCluster(numberOfObservations,false);
  int currentClusterLabel = uniqClusterLabels[clusterNum];

  //pick out observations in the current cluster.
  //clusteringIn.size() == dataValsIn[0].size().
  //that is, the clustering equals the number of rows in the initial data set.
  for (auto currentObs = 0; currentObs != clusteringIn.size(); ++currentObs) { 
    if (clusteringIn[currentObs] == currentClusterLabel)  {
      inCluster[currentObs] = true;
    }
  }
  
  //using the subset flags, get the correspond data from dataValsIn
  std::vector<std::vector<double>> dvSubset;  
  dvSubset.resize(numberOfColumns);
  std::vector<unsigned long> clusterLookup;
    
  for (auto currentObsIndex = 0; currentObsIndex != numberOfObservations; ++currentObsIndex)  {
    if (inCluster[currentObsIndex]) {
      clusterLookup.push_back(currentObsIndex);
      for (auto currentColumn = 0; currentColumn != numberOfColumns; ++currentColumn) {
	(dvSubset[currentColumn]).push_back((dataValsIn[currentColumn])[currentObsIndex]);
      }
    }
  }
  std::vector<double> clusterLowerQs, clusterMedians, clusterUpperQs;
  std::vector<double> qnt = {0.0,0.0,0.0};
  std::vector<double> quantileProbs = annComparisonQuantiles;

  for (auto currentColumn = 0; currentColumn != numberOfColumns; ++currentColumn) {
    qnt = rQuantile(dvSubset[currentColumn],quantileProbs);
    clusterLowerQs.push_back(qnt[0]);
    clusterMedians.push_back(qnt[1]);
    clusterUpperQs.push_back(qnt[2]);
  }
  std::vector<int> newLabel;
  double cqLower,cq50,cqUpper,cDataValue,cAnnotationValue;
  std::vector<double> colAnnBdry;
  int colAnnScore = 0;
  unsigned int numAnnPts;

  for (auto clusterObsIndex = 0; clusterObsIndex != clusterLookup.size(); ++clusterObsIndex) {
    newLabel.clear();
    for (auto currentColNum = 0; currentColNum != numberOfColumns; ++currentColNum) {
      if (annotatingCols[currentColNum]) {
	cqLower = clusterLowerQs[currentColNum];
	cq50 = clusterMedians[currentColNum];
	cqUpper = clusterUpperQs[currentColNum];
	colAnnBdry = annotationBoundaries[currentColNum];
	numAnnPts = colAnnBdry.size();
	colAnnScore = 0;
	if (numAnnPts == 1) {
	  //if there is only one annotation boundary, 
	  //the decision to split the cluster follows by compare the lower and upper cluster quantiles to
	  //the the columns annotation boundary. 
	  cDataValue = (dataValsIn[currentColNum])[(clusterLookup[clusterObsIndex])];
	  cAnnotationValue = colAnnBdry[0]; //isCustomLabel == true implies this is the unique annotation value.
	  if (cqLower >= cAnnotationValue) {
	    //if the lower quantile exceeds the single annotation boundary, all obs in cluster == high
	    colAnnScore = 1;
	  }
	  else if (cqUpper <= cAnnotationValue) {
	    //if the upper quantile falls below the single annotation boundary, all obs in cluster == low
	    colAnnScore = 0;
	  }
	  else if (cDataValue >= cAnnotationValue) {
	    //if both previous conditions fail, split the cluster.
	    colAnnScore = 1;
	  }
	  else {
	    colAnnScore = 0;
	  }
	}
	else {
	  //with multiple annotation boundaries, simply use the median.
	  for (auto j = 0; j != numAnnPts; ++j) {
	    if (cq50 >= colAnnBdry[j]){
	      colAnnScore += 1;
	    }
	  }
	}
	//derive annotations for output.
	newLabel.push_back(colAnnScore);
	newLabel.push_back(numAnnPts);
      }
      else {
	newLabel.push_back(0);
	newLabel.push_back(0);
      }
    }
    //update annotations with the new label
    annotatedClusters[(clusterLookup[clusterObsIndex])] = newLabel;
  }

  //lock the mutux to update the annotation store
  std::lock_guard<std::mutex> updateAnn(annMutRef);
  std::vector<int> parUpdate;
  for (auto nlv = 0; nlv != annotatedClusters.size(); ++nlv) {
    parUpdate = annotatedClusters[nlv];
    //member of a cluster will be twice the numberOfColumns.
    //all other entries in a thread will be exactly the numberOfColumns,
    //which we do not need to copy to the global annotated vector.
    if (parUpdate.size() != numberOfColumns) {
      globalAnnotatedClusters[nlv] = parUpdate;
    }
  }
  return;
}


std::vector<std::vector<int>> parallelAnnotateCluster(const std::vector<std::vector<double>>& dataValsIn,
						      const std::vector<int>& clusteringIn,
						      const std::vector<bool>& annotatingCols,
						      const std::vector<std::vector<double>>& annotationBoundaries,
						      const std::vector<std::string>& annotationsIn,
						      const std::vector<double>& annComparisonQuantiles,
						      const unsigned long numThreadsInUse) { 

  //first determine the unique cluster labels in the input clustering
  std::vector<int> uniqClusterLabels = clusteringIn;
  std::sort(uniqClusterLabels.begin(),uniqClusterLabels.end());
  auto uit = std::unique(uniqClusterLabels.begin(),uniqClusterLabels.end());
  uniqClusterLabels.resize(std::distance(uniqClusterLabels.begin(),uit));
  unsigned long numberOfClusters = uniqClusterLabels.size();

  //make a template to store annotated clusters.
  std::vector<int> labelTemplate((2*(dataValsIn.size())),0);
  std::vector<std::vector<int>> globalAnnotatedClusters((clusteringIn.size()),labelTemplate);

  unsigned long clusterNum = 0;
  std::mutex annMutex;
  std::vector<std::thread> threads;
  while (clusterNum != numberOfClusters) {
    threads.clear();  //remove existing threads (if any)
    for (auto i = 0; i < (numThreadsInUse-1); ++i) {
      if (clusterNum < numberOfClusters) {
	//if there are still clusters to annotate, launch a thread to annotate it.
	threads.push_back(std::thread(threadAnnotateCluster,
				      std::ref(dataValsIn),
				      std::ref(clusteringIn),
				      std::ref(annotatingCols),
				      std::ref(annotationBoundaries),
				      std::ref(annotationsIn),
				      std::ref(annComparisonQuantiles),
				      std::ref(uniqClusterLabels),
				      clusterNum,
				      std::ref(annMutex),
				      std::ref(globalAnnotatedClusters)));
	++clusterNum;
      }
      else {
	//we have annotated every cluster
	continue;
      }
    }
    //join all threads we launched
    std::for_each(threads.begin(),threads.end(),
		  std::mem_fn(&std::thread::join));
  }
  return globalAnnotatedClusters;
}
