.constructAnalysisMap <- function(
                                  projectPath,
                                  gspData,
                                  sampNames,
                                  experimentalUnit,
                                  imputationHierarchy,
                                  debugFlag
                                  )
{
    #construct the analysis map directly from gating set pData
    if ((experimentalUnit == "") || (!(experimentalUnit %in% colnames(gspData)))) {
        analysisMap <- data.frame(
            sampleName = sampNames,
            experimentalUnit = sampNames,
            stringsAsFactors = FALSE
        )
    }
    else {
        analysisMap <- data.frame(
            sampleName = sampNames,
            experimentalUnit = gspData[,experimentalUnit,drop=TRUE],
            stringsAsFactors = FALSE
        )
    }
    if ((imputationHierarchy=="") || (!(imputationHierarchy %in% colnames(gspData)))) {
        analysisMap$impH <- "allSamples"
    }
    else {
        analysisMap$impH <- as.character(gspData[,imputationHierarchy,drop=TRUE])
    }
    if (debugFlag) {
        print("impH")
        print(table(analysisMap$impH))
    }
    #test to see if the analysis map has changed between FAUST runs
    if (!file.exists(file.path(normalizePath(projectPath),
                               "faustData",
                               "metaData",
                               "analysisMap.rds"))) {
        saveRDS(analysisMap,file.path(normalizePath(projectPath),
                                      "faustData",
                                      "metaData",
                                      "analysisMap.rds"))
    }
    else {
        oldAnalysisMap <- readRDS(file.path(normalizePath(projectPath),
                                            "faustData",
                                            "metaData",
                                            "analysisMap.rds"))
        if (!identical(oldAnalysisMap,analysisMap)) {
            if (nrow(oldAnalysisMap) != nrow(analysisMap)) {
                print("The number of samples has changed between faust runs.")
                stop("Please start a new projectPath to analyze a different collection of samples.")
            }
            else {
                print("The sample grouping derived from experimentalUnit has changed between faust runs.")
                stop("Please start a new projectPath to analyze a different concatenation of samples.")
            }
        }
    }
    return()
}
