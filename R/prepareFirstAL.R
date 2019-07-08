.prepareFirstAL <- function(analysisMap,projectPath) {
    if (!dir.exists(paste0(projectPath,"/faustData/levelData"))) {
        dir.create(paste0(projectPath,"/faustData/levelData"))
    }
    if (!file.exists(paste0(projectPath,"/faustData/metaData/firstALReady.rds"))) {
        uniqueLevels <- unique(analysisMap[,"analysisLevel"])
        for (aLevel in uniqueLevels) {
            aData <- analysisMap[which(analysisMap[,"analysisLevel"]==aLevel),,drop=FALSE]
            firstSample <- TRUE
            for (sampleNum in seq(nrow(aData))) {
                sampleName <- aData[sampleNum,"sampleName"]
                if (firstSample) {
                    levelExprs <- readRDS(paste0(projectPath,"/faustData/sampleData/",sampleName,"/exprsMat.rds"))
                    levelRes <- readRDS(paste0(projectPath,"/faustData/sampleData/",sampleName,"/resMat.rds"))
                    levelLookup <- rep(sampleName,nrow(levelExprs))
                    firstSample <- FALSE
                }
                else {
                    newExprs <- readRDS(paste0(projectPath,"/faustData/sampleData/",sampleName,"/exprsMat.rds"))
                    levelExprs <- rbind(levelExprs,newExprs)
                    newRes <- readRDS(paste0(projectPath,"/faustData/sampleData/",sampleName,"/resMat.rds"))
                    levelRes <- rbind(levelRes,newRes)
                    newLookup <- rep(sampleName,nrow(newExprs))
                    levelLookup <- append(levelLookup,newLookup)
                }
            }
            if (nrow(levelExprs)) { #there is data associated with the analysis level. record it.
                if (!dir.exists(paste0(projectPath,"/faustData/levelData/",aLevel))) {
                    dir.create(paste0(projectPath,"/faustData/levelData/",aLevel))
                }
                saveRDS(levelExprs,paste0(projectPath,"/faustData/levelData/",aLevel,"/levelExprs.rds"))
                saveRDS(levelRes,paste0(projectPath,"/faustData/levelData/",aLevel,"/levelRes.rds"))
                saveRDS(levelLookup,paste0(projectPath,"/faustData/levelData/",aLevel,"/levelLookup.rds"))
            }
        }
        firstALReady <- TRUE
        saveRDS(firstALReady,paste0(projectPath,"/faustData/metaData/firstALReady.rds"))
    }
    return()
}
