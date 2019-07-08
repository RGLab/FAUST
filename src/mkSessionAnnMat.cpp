#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <type_traits>
#include <iomanip>
#include <Rcpp.h>

template <typename T>
std::vector<T> mySplitStr(std::string inString, std::string delim, T typeSignal)
{
  std::string str2parse = inString;
  std::vector<T> outVec;
  size_t pos = 0;
  std::string subStr;
  while ((pos = str2parse.find(delim)) != std::string::npos) {
    subStr = str2parse.substr(0, pos);
    if (std::is_same<T,int>::value)
      outVec.push_back(std::stoi(subStr));
    else
      outVec.push_back(std::stod(subStr));
    str2parse.erase(0, pos + delim.length());
  }
  //push back remaining element.
  if (std::is_same<T,int>::value)
    outVec.push_back(std::stoi(str2parse));
  else
    outVec.push_back(std::stod(str2parse));
  return outVec;
}


void updateAnnMat(std::string expMatFP,
		  std::string newAnnMatFP,
		  std::vector<std::vector<double>> cutPoints) {
  std::string expLine, newAnn;
  //std::ifstream annFile(annMatFP,std::ifstream::in);
  std::ifstream expFile(expMatFP,std::ifstream::in);
  std::ofstream newAnnFile(newAnnMatFP);
  std::vector<int> av; 
  int annotationScore;
  std::vector<double> ev, cp;
  while (std::getline(expFile,expLine)) {
    ev = mySplitStr(expLine,",",2.01);
    newAnn = "";
    for (auto colNum = 0; colNum != cutPoints.size();++colNum) {
      cp = cutPoints[colNum];
      annotationScore = 1;
      for (auto cpNum =0; cpNum != cp.size(); ++cpNum) {
	if (ev[colNum] >= cp[cpNum]) {
	  ++annotationScore;
	}
      }
      newAnn += std::to_string(annotationScore);
      if (colNum < (cutPoints.size() -1))
	newAnn += ",";
    }
    newAnn += "\n";
    newAnnFile << newAnn;
  }
  expFile.close();
  newAnnFile.close();
  return;
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
void mkSessionAnnMat(std::string expMatFP,
		     std::string newAMFP,
		     Rcpp::List cpList) {
  //copy cut points into cpp data structure
  std::vector<std::vector<double>> cutPoints;
  std::vector<double> tmpVec;
  for (auto ncp = 0; ncp != cpList.size(); ++ncp) {
    tmpVec = Rcpp::as<std::vector<double>>(cpList[ncp]);
    cutPoints.push_back(tmpVec);
  }
  //write a new annotation matrix to file.
  updateAnnMat(expMatFP,newAMFP,cutPoints);
  return;
}
