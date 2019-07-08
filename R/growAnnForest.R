.growForestForALevel <- function(aLevel,
                                 rootPop,
                                 activeChannels,
                                 analysisMap,
                                 numIter,
                                 debugFlag,
                                 threadNum,
                                 seedValue,
                                 projectPath)
{
    #function to 
    if (debugFlag) print(paste0("Growing annotation forest for: ",aLevel))
    levelExprs <- readRDS(paste0(projectPath,"/faustData/levelData/",aLevel,"/levelExprs.rds"))
    levelRes <- readRDS(paste0(projectPath,"/faustData/levelData/",aLevel,"/levelRes.rds"))
    levelExprs <- levelExprs[,activeChannels,drop=FALSE]
    levelRes <- levelRes[,activeChannels,drop=FALSE]
    resFlag <- FALSE
    for (colNum in seq(ncol(levelRes))) {
        if (length(which(levelRes[,colNum] > 0))) {
            resFlag <- TRUE
            break
        }
    }
    annF <- growAnnotationForest(
        dataSet=levelExprs,
        numberIterations=numIter,
        pValueThreshold=0.25,
        minimumClusterSize=25,
        randomCandidateSearch=FALSE,
        maximumSearchDepth=2,
        numberOfThreads=threadNum,
        maximumGatingNum=1e10,
        anyValueRestricted=resFlag,
        resValMatrix=levelRes,
        cutPointUpperBound=2,
        getDebugInfo=FALSE,
        randomSeed=seedValue,
        recordCounts=FALSE,
        recordIndices=FALSE
    )
    saveRDS(annF,paste0(projectPath,"/faustData/levelData/",aLevel,"/",rootPop,"_annF.rds"))
    ePop <- apply(levelRes,2,function(x){length(which(x==0))})
    names(ePop) <- colnames(levelExprs)
    af <- annF[["gateData"]]
    pAnnF <- .parseAnnotationForest(af,ePop)
    saveRDS(pAnnF,paste0(projectPath,"/faustData/levelData/",aLevel,"/",rootPop,"_pAnnF.rds"))
    aLevelDone <- TRUE
    saveRDS(aLevelDone,paste0(projectPath,"/faustData/levelData/",aLevel,"/aLevelComplete.rds"))
    if (debugFlag) print(paste0("Annotation forest complete for: ",aLevel))
    return()
}



.growAnnForest <- function(rootPop,
                           activeChannels,
                           analysisMap,
                           numIter,
                           debugFlag,
                           threadNum,
                           seedValue,
                           projectPath)
{
    uniqueLevels <- sort(unique(analysisMap[,"analysisLevel"]))
    activeLevels <- c()
    #accumulate vector of levels without annotation forests.
    for (analysisLevel in uniqueLevels) {
        if (!file.exists(paste0(projectPath,"/faustData/levelData/",analysisLevel,"/aLevelComplete.rds"))) {
            activeLevels <- append(activeLevels,analysisLevel)
        }
    }
    #grow forests for levels that lack them
    while (length(activeLevels)) {
        currentLevel <- activeLevels[1]
        activeLevels <- activeLevels[-1]
        .growForestForALevel(
            aLevel=currentLevel,
            rootPop=rootPop,
            activeChannels=activeChannels,
            analysisMap=analysisMap,
            numIter=numIter,
            debugFlag=debugFlag,
            threadNum=threadNum,
            seedValue=seedValue,
            projectPath=projectPath
        )
    }
    return()
}

