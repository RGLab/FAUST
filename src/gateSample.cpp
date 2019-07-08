#include <Rcpp.h>
#include <string>
#include <sstream>
#include <unordered_set>
// An Rcpp gateSample function. 
using namespace Rcpp;
// [[Rcpp::plugins("cpp11")]]

// [[Rcpp::export]]
StringVector gateSample(NumericMatrix annotationMatrix, StringVector selectedChannels, NumericVector gateNums, StringVector scampCellPops) {
  int nr = annotationMatrix.nrow();
  int nc = annotationMatrix.ncol();
  std::ostringstream oss;
  Rcpp::StringVector resultvector(nr); // the returned vector.
  std::unordered_set<std::string> allowed_pops;
  // add the allowed annotations to an unordered set.
  for(auto& s : scampCellPops){
    allowed_pops.insert(Rcpp::as< std::string >(s));
  }
  //this insertion doesn't compile on the rhinos
  // allowed_pops.insert(scampCellPops.begin(), scampCellPops.end());
  // loop through each annotation for each cell and construct the phenotype string.
  for(int r = 0; r < nr; r++){
    for(int c = 0; c < nc; c++){
      oss << selectedChannels(c) << "~" << annotationMatrix(r,c) << "~" << gateNums(c) << "~";
    }
    std::string v = oss.str();
   
    // test if it's a valid annotation against the unordered set.
    std::unordered_set<std::string>::const_iterator got = allowed_pops.find(v);
    if(got != allowed_pops.end()){ // currently not finding the string.
      // found the string so add the annotation to the results matrix
      resultvector(r) = v;
    }else{
      // didn't find the string. so add 0_0_0_0_0 to the result matrix.
      resultvector(r) = "0_0_0_0_0";
    }
    // clear the output string stream
    oss.str("");
  }
  return resultvector;
}



// This code tests the C++ gating function based on data at the path below.


/*** R
projectPath <- "/shared/silo_researcher/Gottardo_R/gfinak_working/CITN07/FAUST_pheno_full_remapped_persample/"
startingCellPop <- "45"
activeSamples <- list.files(paste0(projectPath,"faustData/sampleData"))
sampleName <- activeSamples[1]
selectedChannels <- readRDS(paste0(projectPath,"faustData/metaData/initSelC.rds"))
analysisMap <- data.frame(sampleName = activeSamples, analysisLevel = activeSamples, stringsAsFactors = FALSE)


{
  resList <- readRDS(paste0(projectPath,"/faustData/gateData/",startingCellPop,"_resList.rds"))
  activeSamples <- analysisMap[,"sampleName"]
  scampCellPops <- readRDS(paste0(projectPath,"/faustData/metaData/scampClusterNames.rds"))
  # for (sampleName in activeSamples) {
  sampleName <- activeSamples[1]
    sAnn <- utils::read.table(file=paste0(projectPath,"/faustData/sampleData/",sampleName,"/scampAnnotation.csv"),
                              header=F,sep="`",
                              stringsAsFactors=FALSE)[,1]
    annotationMatrix <- utils::read.table(file=paste0(projectPath,"/faustData/sampleData/",sampleName,"/annotationMatrix.csv"),
                                          header=F,
                                          sep=",",
                                          stringsAsFactors=FALSE)
    colnames(annotationMatrix) <- selectedChannels
    if (!(length(sAnn)==nrow(annotationMatrix))) {
      print(paste0("Annotation matrix doesn't match scamp annotation in ",sampleName))
      stop("Investigate.")
    }
    aLevel <- analysisMap[which(analysisMap[,"sampleName"]==sampleName),"analysisLevel"]
    scampAF <- rep(list(NA),length(selectedChannels))
    names(scampAF) <- selectedChannels
    for (channel in selectedChannels) {
      scampAF[[channel]] <- resList[[channel]][[aLevel]]
    }
    gateNums <- unlist(lapply(scampAF,length)) + 1
    exactPartition <- rep("0_0_0_0_0",nrow(annotationMatrix))
    c_result <- gateSample(as.matrix(annotationMatrix), selectedChannels, gateNums, scampCellPops);
    for (rowNum in seq(nrow(annotationMatrix))) {
      ep <- paste0(paste0(paste0(selectedChannels,"~",annotationMatrix[rowNum,]),"~",gateNums,"~"),collapse="")
      if (ep %in% scampCellPops) {
        exactPartition[rowNum] <- ep
      }
    }
    all.equal(c_result, exactPartition)
    # data.table::fwrite(list(exactPartition),
    #                    file=paste0(projectPath,"/faustData/sampleData/",sampleName,"/faustAnnotation.csv"),
    #                    sep="~",
    #                    append=FALSE,
    #                    row.names=FALSE,
    #                    col.names=FALSE,
    #                    quote=FALSE)
  # }
  # return()
}
*/
