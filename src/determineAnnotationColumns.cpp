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
#include <math.h> 
#include <algorithm>
#include <numeric>
#include <Rcpp.h>
#include <unordered_map>
#include "faust.h"

std::vector<bool> determineAnnotationColumns(const std::vector<std::vector<double>>& dataValsIn,
					     unsigned long colNum,
					     double dipThreshold,
					     const std::vector<std::vector<int>>&resValsIn,
					     const bool useRestrictedVals){
  std::vector<std::vector<double>> sortedDataVals;
  std::vector<int> restrictObs;
  std::vector<double> valsForP;
  if (useRestrictedVals) {
    // if values are restricted (and so have a score of 1 or 2), omit them from the sorted data.
    // assumption: phenotypical boundaries are determined only by antimodes of **non-zero** cytof observations.
    for (auto i = 0; i != resValsIn.size(); ++i) {
      valsForP.clear();
      restrictObs = resValsIn[i];
      for (auto j = 0; j != restrictObs.size(); ++j) {
	if (restrictObs[j]==0) {
	  valsForP.push_back(((dataValsIn[i])[j]));
	}
      }
      sortedDataVals.push_back(valsForP);
    }
  }
  else {
    //if no values are restricted, use all available data.
    sortedDataVals = dataValsIn;
  }

  //copy & sort the raw data
  for (auto it = sortedDataVals.begin(); it != sortedDataVals.end(); ++it) {
      std::sort((*it).begin(),(*it).end());
  }
  std::vector<bool> labelColumn(colNum,false);
  std::vector<double> curSortData, depthNum, depthDenom;
  columnSummary curSum;
  for (int i = 0; i != colNum; ++i) {
    curSortData.clear();
    curSortData = sortedDataVals[i];
    if (singleDip(curSortData) < dipThreshold) {
      labelColumn[i]=true;
    }
  }
  return labelColumn;
}
