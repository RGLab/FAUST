.clusterALevelWithScamp <- function(aLevel,
                                    resList,
                                    selectedChannels,
                                    numScampIter,
                                    minClusterSize,
                                    debugFlag,
                                    threadNum,
                                    seedValue,
                                    projectPath
                                    )
{
    if (debugFlag) print(paste0("Starting SCAMP for: ",aLevel))
    levelExprs <- readRDS(paste0(projectPath,"/faustData/levelData/",aLevel,"/levelExprs.rds"))
    levelRes <- readRDS(paste0(projectPath,"/faustData/levelData/",aLevel,"/levelRes.rds"))
    levelExprs <- levelExprs[,selectedChannels, drop = FALSE]
    levelRes <- levelRes[,selectedChannels, drop = FALSE]
    resFlag <- FALSE
    for (colNum in seq(ncol(levelRes))) {
        if (length(which(levelRes[,colNum] > 0))) {
            resFlag <- TRUE
            break
        }
    }
    #get the scamp annotation forest
    scampAF <- rep(list(NA),ncol(levelExprs))
    names(scampAF) <- colnames(levelExprs)
    for (colName in colnames(levelExprs)) {
        scampAF[[colName]] <- resList[[colName]][[aLevel]]
    }
    #scale the minimum cluster size to respect the number of sample in the analysisLevel.
    if (minClusterSize > 0) {
        useClusterNum <- minClusterSize
    }
    else {
        useClusterNum <- 100 * length(uniqueLevels)
    }
    scampClustering <- scamp(
        dataSet = levelExprs,
        numberIterations = numScampIter,
        minimumClusterSize = useClusterNum,
        numberOfThreads = threadNum,
        anyValueRestricted = resFlag,
        resValMatrix = levelRes,
        useAnnForest = TRUE,
        annForestVals = scampAF,
        randomSeed=seedValue,
        getDebugInfo = FALSE
    )
    runClustering <- scampClustering[[1]]
    maxClustering <- scampClustering[[2]]
    outClustering <- rep("Uncertain",length(maxClustering))
    agreeIndex <- which(runClustering==maxClustering)
    if (length(agreeIndex)) {
        outClustering[agreeIndex] <- maxClustering[agreeIndex]
    }
    else {
        print(paste0("SCAMP reached no concensus for: ",aLevel))
        print("This indicates high-amount of uncertainty in cluster labels.")
        print("Increase the numScampIter parameters until concensus reached")
        stop("killing job.")
    }
    clusterNames <- setdiff(sort(unique(names(table(outClustering)))),"Uncertain")
    saveRDS(clusterNames,paste0(projectPath,"/faustData/levelData/",aLevel,"/scampClusterLabels.rds"))
    #unwind the level to each sample
    levelLookup <- readRDS(paste0(projectPath,"/faustData/levelData/",aLevel,"/levelLookup.rds"))
    for (sampleName in names(table(levelLookup))) {
        sampleLookup <- which(levelLookup == sampleName)
        if (length(sampleLookup)) {
            sampleClustering <- outClustering[sampleLookup]
            data.table::fwrite(list(sampleClustering),
                               file = paste0(projectPath,"/faustData/sampleData/",sampleName,"/scampAnnotation.csv"),
                               sep = "`",
                               append = FALSE,
                               row.names = FALSE,
                               col.names = FALSE,
                               quote = FALSE)
        }
    }
    scampALevelDone <- TRUE
    saveRDS(scampALevelDone,paste0(projectPath,"/faustData/levelData/",aLevel,"/scampALevelComplete.rds"))
    if (debugFlag) print(paste0("SCAMP complete for: ",aLevel))
    return()
}

                                    
.clusterLevelsWithScamp <- function(startingCellPop,
                                    selectedChannels,
                                    analysisMap,
                                    numScampIter,
                                    nameOccuranceNum,
                                    minClusterSize,
                                    debugFlag,
                                    threadNum,
                                    seedValue,
                                    projectPath
                                    )
{
    resList <- readRDS(paste0(projectPath,"/faustData/gateData/",startingCellPop,"_resList.rds"))
    uniqueLevels <- sort(unique(analysisMap[,"analysisLevel"]))
    activeLevels <- c()
    #accumulate vector of levels without annotation forests.
    for (analysisLevel in uniqueLevels) {
        if (!file.exists(paste0(projectPath,"/faustData/levelData/",analysisLevel,"/scampALevelComplete.rds"))) {
            activeLevels <- append(activeLevels,analysisLevel)
        }
    }
    #grow forests for levels that lack them
    while (length(activeLevels)) {
        currentLevel <- activeLevels[1]
        activeLevels <- activeLevels[-1]
        .clusterALevelWithScamp(
            aLevel=currentLevel,
            resList=resList,
            selectedChannels=selectedChannels,
            numScampIter=numScampIter,
            minClusterSize=minClusterSize,
            debugFlag=debugFlag,
            threadNum=threadNum,
            seedValue=seedValue,
            projectPath=projectPath
        )
    }
    #collect all labels from the scamp clusterings
    if (debugFlag) print("Accumulating cluster labels.")
    clusterNames <- c()
    for (analysisLevel in uniqueLevels) {
        if (!file.exists(paste0(projectPath,"/faustData/levelData/",analysisLevel,"/scampClusterLabels.rds"))) {
            print(paste0("Labels not detected in analysisLevel ",analysisLevel))
            print("This is a bug -- all analysis levels should have labels.")
            stop("Killing FAUST. Check logs to determine which level is unlabeled.")
        }
        else {
            levelLabels <- readRDS(paste0(projectPath,"/faustData/levelData/",analysisLevel,"/scampClusterLabels.rds"))
            clusterNames <- append(clusterNames,levelLabels)
        }
    }
    nameSummary <- table(clusterNames)
    saveRDS(nameSummary,paste0(projectPath,"/faustData/metaData/scampNameSummary.rds"))
    clusterNames <- names(nameSummary[which(nameSummary >= nameOccuranceNum)])
    saveRDS(clusterNames,paste0(projectPath,"/faustData/metaData/scampClusterNames.rds"))
    nameSummaryPlotDF <- data.frame(x=seq(max(nameSummary)),
                                    y=sapply(seq(max(nameSummary)),function(x){
                                        length(which(nameSummary >= x))}))
    nspOut <- ggplot(nameSummaryPlotDF,aes(x=x,y=y))+
        geom_line()+
        theme_bw()+
        geom_vline(xintercept=nameOccuranceNum,col="red")+
        xlab("Number of times a cluster name appears across SCAMP clusterings")+
        ylab("Number of SCAMP clusters >= the appearance number")+
        ggtitle("Red line is nameOccuranceNum setting in faust")
    cowplot::save_plot(paste0(projectPath,"/faustData/plotData/scampNamesPlot.pdf"),
                       nspOut,base_height=15,base_width=15)
    if (debugFlag) print("Cluster labels collected and saved.")
    return()
}
