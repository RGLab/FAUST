.prepareFirstAL <- function(projectPath) {
    if (!dir.exists(file.path(normalizePath(projectPath),
                              "faustData",
                              "levelData")))
    {
        dir.create(file.path(normalizePath(projectPath),
                             "faustData",
                             "levelData"))
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
        uniqueLevels <- unique(analysisMap[,"analysisLevel"])
        for (aLevel in uniqueLevels) {
            aData <- analysisMap[which(analysisMap[,"analysisLevel"]==aLevel),,drop=FALSE]
            firstSample <- TRUE
            for (sampleNum in seq(nrow(aData))) {
                sampleName <- aData[sampleNum,"sampleName"]
                if (firstSample) {
                    levelExprs <- readRDS(file.path(normalizePath(projectPath),
                                                    "faustData",
                                                    "sampleData",
                                                    sampleName,
                                                    "exprsMat.rds"))
                    levelRes <- readRDS(file.path(normalizePath(projectPath),
                                                  "faustData",
                                                  "sampleData",
                                                  sampleName,
                                                  "resMat.rds"))
                    levelLookup <- rep(sampleName,nrow(levelExprs))
                    firstSample <- FALSE
                }
                else {
                    newExprs <- readRDS(file.path(normalizePath(projectPath),
                                                  "faustData",
                                                  "sampleData",
                                                  sampleName,
                                                  "exprsMat.rds"))
                    levelExprs <- rbind(levelExprs,newExprs)
                    newRes <- readRDS(file.path(normalizePath(projectPath),
                                                "faustData",
                                                "sampleData",
                                                sampleName,
                                                "resMat.rds"))
                    levelRes <- rbind(levelRes,newRes)
                    newLookup <- rep(sampleName,nrow(newExprs))
                    levelLookup <- append(levelLookup,newLookup)
                }
            }
            if (nrow(levelExprs)) { #there is data associated with the analysis level. record it.
                if (!dir.exists(file.path(normalizePath(projectPath),
                                          "faustData",
                                          "levelData",
                                          aLevel))) {
                    dir.create(file.path(normalizePath(projectPath),
                                         "faustData",
                                         "levelData",
                                         aLevel))
                }
                saveRDS(levelExprs,
                        file.path(normalizePath(projectPath),
                                  "faustData",
                                  "levelData",
                                  aLevel,
                                  "levelExprs.rds"))
                saveRDS(levelRes,
                        file.path(normalizePath(projectPath),
                                  "faustData",
                                  "levelData",
                                  aLevel,
                                  "levelRes.rds"))
                saveRDS(levelLookup,
                        file.path(normalizePath(projectPath),
                                  "faustData",
                                  "levelData",
                                  aLevel,
                                  "levelLookup.rds"))
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
