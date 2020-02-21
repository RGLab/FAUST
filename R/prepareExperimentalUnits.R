.prepareExperimentalUnits <- function(projectPath) {
    if (!dir.exists(file.path(normalizePath(projectPath),
                              "faustData",
                              "expUnitData")))
    {
        dir.create(file.path(normalizePath(projectPath),
                             "faustData",
                             "expUnitData"))
    }
    if (!file.exists(file.path(normalizePath(projectPath),
                               "faustData",
                               "metaData",
                               "firstALReady.rds")))
    {
        analysisMap <- readRDS(file.path(normalizePath(projectPath),
                                         "faustData",
                                         "metaData",
                                         "analysisMap.rds"))
        uniqueExpUnits <- unique(analysisMap[,"experimentalUnit"])
        for (expUnit in uniqueExpUnits) {
            aData <- analysisMap[which(analysisMap[,"experimentalUnit"]==expUnit),,drop=FALSE]
            firstSample <- TRUE
            for (sampleNum in seq(nrow(aData))) {
                sampleName <- aData[sampleNum,"sampleName"]
                if (firstSample) {
                    expUnitExprs <- readRDS(file.path(normalizePath(projectPath),
                                                    "faustData",
                                                    "sampleData",
                                                    sampleName,
                                                    "exprsMat.rds"))
                    expUnitRes <- readRDS(file.path(normalizePath(projectPath),
                                                  "faustData",
                                                  "sampleData",
                                                  sampleName,
                                                  "resMat.rds"))
                    expUnitToSampleLookup <- rep(sampleName,nrow(expUnitExprs))
                    firstSample <- FALSE
                }
                else {
                    newExprs <- readRDS(file.path(normalizePath(projectPath),
                                                  "faustData",
                                                  "sampleData",
                                                  sampleName,
                                                  "exprsMat.rds"))
                    expUnitExprs <- rbind(expUnitExprs,newExprs)
                    newRes <- readRDS(file.path(normalizePath(projectPath),
                                                "faustData",
                                                "sampleData",
                                                sampleName,
                                                "resMat.rds"))
                    expUnitRes <- rbind(expUnitRes,newRes)
                    newLookup <- rep(sampleName,nrow(newExprs))
                    expUnitToSampleLookup <- append(expUnitToSampleLookup,newLookup)
                }
            }
            if (nrow(expUnitExprs)) { #there is data associated with the experimental unit. record it.
                if (!dir.exists(file.path(normalizePath(projectPath),
                                          "faustData",
                                          "expUnitData",
                                          expUnit))) {
                    dir.create(file.path(normalizePath(projectPath),
                                         "faustData",
                                         "expUnitData",
                                         expUnit))
                }
                saveRDS(expUnitExprs,
                        file.path(normalizePath(projectPath),
                                  "faustData",
                                  "expUnitData",
                                  expUnit,
                                  "expUnitExprs.rds"))
                saveRDS(expUnitRes,
                        file.path(normalizePath(projectPath),
                                  "faustData",
                                  "expUnitData",
                                  expUnit,
                                  "expUnitRes.rds"))
                saveRDS(expUnitToSampleLookup,
                        file.path(normalizePath(projectPath),
                                  "faustData",
                                  "expUnitData",
                                  expUnit,
                                  "expUnitToSampleLookup.rds"))
            }
        }
        firstALReady <- TRUE
        saveRDS(firstALReady,
                file.path(normalizePath(projectPath),
                          "faustData",
                          "metaData",
                          "firstALReady.rds"))
    }
    return()
}
