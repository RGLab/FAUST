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

#include <string>
#include <vector>
#include <Rcpp.h>
#include <iostream>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
std::vector<std::string> labelPartition(const Rcpp::NumericMatrix& partitionScores,
					const std::vector<std::string>& colNames,
					const std::vector<int>& colScoreRange)
{
  //default label containers.
  std::vector<std::vector<std::string>> finalLabels(9);
  std::vector<std::string> oneLabel = {"NoCut"};
  std::vector<std::string> twoLabels = {"Lowest","Highest"};
  std::vector<std::string> threeLabels = {"Lowest","Medium","Highest"};
  std::vector<std::string> fourLabels = {"Lowest","MediumLow","MediumHigh","Highest"};
  std::vector<std::string> fiveLabels = {"Lowest","MediumLow","Medium","MediumHigh","Highest"};
  std::vector<std::string> sixLabels = {"Lowest","Low","MediumLow","MediumHigh","High","Highest"};
  std::vector<std::string> sevenLabels = {"Lowest","Low","MediumLow","Medium",
					  "MediumHigh","High","Highest"};
  std::vector<std::string> eightLabels = {"Lowest","VeryLow","Low","MediumLow","MediumHigh",
					  "High","VeryHigh","Highest"};
  std::vector<std::string> nineLabels = {"Lowest","VeryLow","Low","MediumLow","Medium",
					 "MediumHigh","High","VeryHigh","Highest"};
  finalLabels[0]=oneLabel;
  finalLabels[1]=twoLabels;
  finalLabels[2]=threeLabels;
  finalLabels[3]=fourLabels;
  finalLabels[4]=fiveLabels;
  finalLabels[5]=sixLabels;
  finalLabels[6]=sevenLabels;
  finalLabels[7]=eightLabels;
  finalLabels[8]=nineLabels;

  //construct labels
  std::string annotationString = "";
  std::string tmpStr;
  auto nrows =  partitionScores.nrow();
  auto ncols =  partitionScores.ncol();
  std::vector<std::string> labelVec(nrows,annotationString);
  int adjPscore = 0;
  for (auto i=0; i != nrows; ++i) {
    //Rcpp::Rcout << i << std::endl;
    annotationString = "";
    for (auto j=0; j != ncols; ++j) {
      adjPscore = (partitionScores(i,j)-1);
      if (adjPscore < 0) 
	adjPscore = 0; //the column is labeled NA
      if (j < (ncols -1))
	tmpStr = colNames[j] + "_" + (finalLabels[colScoreRange[j]])[adjPscore] + "_";
      else
	tmpStr = colNames[j] + "_" + (finalLabels[colScoreRange[j]])[adjPscore]; 
      annotationString += tmpStr;
    }
    labelVec[i]= annotationString;
  }
  return labelVec;
}
