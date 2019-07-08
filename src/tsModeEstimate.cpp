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
#include "faust.h"
#include "kMedDP.h"


bool isLocalMaxTS(double &x_cand, double &x_left, double &x_right) {
  //inequality strict because the string jumps to discrete values
  if ((x_cand > x_left) && (x_cand > x_right))  return true;
  else return false;
}


std::vector<double> tsModeEstimate(const std::vector<double> &xVec) {
  //Given sorted univarite data vector xVec, returns an estimate of cut points to classify data
  //Assumes the DIP test has failed to reject the null hypothesis of unimodality.
  //Consequently, the cut point selected is the location of the mode itself.
    
  //In the event that the taut string estimates >1 mode while the DIP rejects the hypothesis of unimodality,
  //return the sample median as the cut point.
  stringInfo tautString = cpPmden(xVec); 
  int numModes = tautString.nmax-(tautString.nmax/2); //always check to see if the taut string sees unimodality.
                                                      //Note we are deliberately taking advantage of integer truncation
                                                      //No need to call std::floor.
  std::vector<double> qProbVal = {0.5};
  std::vector<double> sampQnt;
  double sampleMedian;
  sampQnt = rQuantile(xVec,qProbVal);
  sampleMedian = sampQnt[0];

  if (numModes > 1) {
    Rcpp::Rcout << "Taut string disagrees with DIP; returning sample median." << std::endl;
    std::vector<double> modeGates;
    modeGates.push_back(*std::min_element(xVec.begin(),xVec.end()));
    modeGates.push_back(sampleMedian);
    modeGates.push_back(*std::max_element(xVec.begin(),xVec.end()));
    std::sort(modeGates.begin(),modeGates.end()); //ensure modeGates are returned in ascending order
    return modeGates;
  }
  else {
    std::vector<double> ys = tautString.string;
    std::vector<double> yvals = ys; //copy string into yvals
    auto it = std::unique(yvals.begin(),yvals.end()); //collapse to unique elements
    yvals.resize(std::distance(yvals.begin(),it)); //resize vector
    std::vector<bool> localMax(yvals.size(),false);
    for (auto i = 1; i != (localMax.size()-1); ++i) {
      localMax[i] = isLocalMaxTS(yvals[i],yvals[(i-1)],yvals[(i+1)]);
    }
    
    std::vector<double> cutValue;
    int cutCounter = 0;
    for (auto i = 0; i != yvals.size(); ++i) 
      if (localMax[i] == true) {
	cutValue.push_back(yvals[i]);
	cutCounter += 1;
      }

    std::vector<double> modeGates;
    std::vector<decltype(ys.size())> yInds;
    decltype(ys.size()) cutIndex;
    long long indexSum = 0;
    double cCut; 

    if (cutCounter == 0) {
      Rcpp::Rcout << "Warning: tsModeEstimate string found no mode. Returning sample median." << std::endl;
      modeGates.push_back(*std::min_element(xVec.begin(),xVec.end()));
      modeGates.push_back(sampleMedian);
      modeGates.push_back(*std::max_element(xVec.begin(),xVec.end()));
      std::sort(modeGates.begin(),modeGates.end()); //ensure modeGates are returned in ascending order
      return modeGates;

    }
    else {
      if (cutCounter > 1) {
	Rcpp::Rcout << "Warning: tsModeEstimate finding multiple modes. Using first as cut-point." << std::endl;
      }
      yInds.clear();
      cCut = cutValue[0];
      for (auto i = 0; i != ys.size(); ++i) 
	if (ys[i] == cCut)
	  yInds.push_back(i);
      indexSum = 0;
      for (auto indexVal : yInds)
	indexSum += indexVal;
      cutIndex = indexSum/yInds.size();
    
      modeGates.push_back(*std::min_element(xVec.begin(),xVec.end()));
      modeGates.push_back(xVec[cutIndex]);
      modeGates.push_back(*std::max_element(xVec.begin(),xVec.end()));
      std::sort(modeGates.begin(),modeGates.end()); //ensure modeGates are returned in ascending order
      return modeGates;
    }
  }
}
