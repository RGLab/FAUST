.clusterExpUnitWithScamp <- function(expUnit,
                                    resList,
                                    selectedChannels,
                                    numScampIter,
                                    debugFlag,
                                    threadNum,
                                    seedValue,
                                    projectPath,
                                    densitySubSampleThreshold,
                                    densitySubSampleSize,
                                    densitySubSampleIterations
                                    )
{
    if (debugFlag) print(paste0("Starting SCAMP for: ",expUnit))
    expUnitExprs <- readRDS(file.path(normalizePath(projectPath),
                                    "faustData",
                                    "expUnitData",
                                    expUnit,
                                    "expUnitExprs.rds"))
    expUnitRes <- readRDS(file.path(normalizePath(projectPath),
                                  "faustData",
                                  "expUnitData",
                                  expUnit,
                                  "expUnitRes.rds"))
    expUnitExprs <- expUnitExprs[,selectedChannels, drop = FALSE]
    expUnitRes <- expUnitRes[,selectedChannels, drop = FALSE]
    resFlag <- FALSE
    for (colNum in seq(ncol(expUnitRes))) {
        if (length(which(expUnitRes[,colNum] > 0))) {
            resFlag <- TRUE
            break
        }
    }
    #get the scamp annotation forest
    scampAF <- rep(list(NA),ncol(expUnitExprs))
    names(scampAF) <- colnames(expUnitExprs)
    for (colName in colnames(expUnitExprs)) {
        scampAF[[colName]] <- resList[[colName]][[expUnit]]
    }
    scampClustering <- scamp(
        dataSet = expUnitExprs,
        numberIterations = numScampIter,
        minimumClusterSize = 25,
        numberOfThreads = threadNum,
        anyValueRestricted = resFlag,
        resValMatrix = expUnitRes,
        useAnnForest = TRUE,
        annForestVals = scampAF,
        randomSeed=seedValue,
        getDebugInfo = FALSE,
        subSampleThreshold = densitySubSampleThreshold,
        subSampleSize = densitySubSampleSize,
        subSampleIterations = densitySubSampleIterations
    )
    runClustering <- scampClustering[[1]]
    maxClustering <- scampClustering[[2]]
    outClustering <- rep("Uncertain",length(maxClustering))
    agreeIndex <- which(runClustering==maxClustering)
    if (length(agreeIndex)) {
        outClustering[agreeIndex] <- maxClustering[agreeIndex]
    }
    else {
        print(paste0("SCAMP reached no concensus for: ",expUnit))
        print("This indicates high-amount of uncertainty in cluster labels.")
        print("Increase the numScampIter parameters until concensus reached")
        stop("killing job.")
    }
    clusterNames <- setdiff(sort(unique(names(table(outClustering)))),"Uncertain")
    saveRDS(clusterNames,file.path(normalizePath(projectPath),
                                   "faustData",
                                   "expUnitData",
                                   expUnit,
                                   "scampClusterLabels.rds"))
    #unwind the experimental unit to each sample
    expUnitToSampleLookup <- readRDS(file.path(normalizePath(projectPath),
                                     "faustData",
                                     "expUnitData",
                                     expUnit,
                                     "expUnitToSampleLookup.rds"))
    for (sampleName in names(table(expUnitToSampleLookup))) {
        sampleLookup <- which(expUnitToSampleLookup == sampleName)
        if (length(sampleLookup)) {
            sampleClustering <- outClustering[sampleLookup]
            data.table::fwrite(list(sampleClustering),
                               file = file.path(normalizePath(projectPath),
                                                "faustData",
                                                "sampleData",
                                                sampleName,
                                                "scampAnnotation.csv"),
                               sep = "`",
                               append = FALSE,
                               row.names = FALSE,
                               col.names = FALSE,
                               quote = FALSE)
        }
    }
    scampExpUnitDone <- TRUE
    saveRDS(scampExpUnitDone,
            file.path(normalizePath(projectPath),
                      "faustData",
                      "expUnitData",
                      expUnit,
                      "scampExpUnitComplete.rds"))
    if (debugFlag) print(paste0("SCAMP complete for: ",expUnit))
    return()
}

.computeElbow <- function(phenoDF) {
    if (nrow(phenoDF) < 4) {
        #if a dataset has 3 or fewer samples, default to using all
        #phenotypes from the scamp clusterings.
        elbowLocation <- 1
    }
    else {
        #dynamically estimate the elbow by projecting knots of the GCM 
        #of the distribution of phenotypes onto the line of 
        #observed phenotypes at 5th and 95th rows of the phenoDF
        #trim the extremes to moderate the slopes.
        lowerIndex <- max(2,ceiling(0.05*nrow(phenoDF)))
        upperIndex <- min((nrow(phenoDF)-1),floor(0.95*nrow(phenoDF)))
        elbowDF <- phenoDF[seq(lowerIndex,upperIndex),,drop=FALSE]
        elbX2 <- elbowDF[nrow(elbowDF),"x"]
        elbX1 <- elbowDF[1,"x"]
        elbY2 <- elbowDF[nrow(elbowDF),"y"]
        elbY1 <- elbowDF[1,"y"]
        elbM <- (elbY2-elbY1)/(elbX2-elbX1)
        elbC <- elbY1-(elbM * elbX1)
        lineMin <- function(pt){
            return(abs(pt[2] - elbM*pt[1] - elbC)/sqrt(1+(elbM^2)))
        }
        gcmEst <- fdrtool::gcmlcm(
                               x=elbowDF$x,
                               y=elbowDF$y,
                               type="gcm"
                           )
        candElbow <- apply(cbind(gcmEst$x.knots,gcmEst$y.knots),1,lineMin)
        elbowLocation <- gcmEst$x.knots[which(candElbow==max(candElbow))[1]]
    }
    return(elbowLocation)
}


.clusterExpUnitsWithScamp <- function(projectPath,
                                      nameOccurenceNum,
                                      debugFlag,
                                      threadNum,
                                      seedValue,
                                      densitySubSampleThreshold,
                                      densitySubSampleSize,
                                      densitySubSampleIterations,
                                      archDescriptionList)
{
    #removing numScampIter from interface for simplicity,
    #since it is a rarely modified parameter that makes
    #it harder to apply faust for a new user.
    #if we end up wanting to add numScampIter back to the
    #interface, expose numScampIter = numScampIter
    #in this internal interface.
    numScampIter <- 1

    analysisMap <- readRDS(file.path(normalizePath(projectPath),
                                     "faustData",
                                     "metaData",
                                     "analysisMap.rds"))

    startingCellPop <- readRDS(file.path(normalizePath(projectPath),
                                         "faustData",
                                         "metaData",
                                         "sanitizedCellPopStr.rds"))

    selectedChannels <- readRDS(file.path(normalizePath(projectPath),
                                          "faustData",
                                          "gateData",
                                          paste0(startingCellPop,"_selectedChannels.rds")))

    resList <- readRDS(file.path(normalizePath(projectPath),
                                 "faustData",
                                 "gateData",
                                 paste0(startingCellPop,"_resList.rds")))

    uniqueExpUnits <- sort(unique(analysisMap[,"experimentalUnit"]))
    activeExpUnits <- c()
    #accumulate vector of experimental units without annotation forests.
    for (experimentalUnit in uniqueExpUnits) {
        if (!file.exists(file.path(normalizePath(projectPath),
                                   "faustData",
                                   "expUnitData",
                                   experimentalUnit,
                                   "scampExpUnitComplete.rds"))) {
            activeExpUnits <- append(activeExpUnits,experimentalUnit)
        }
    }
    #grow forests for experimental units that lack them
    if ((length(activeExpUnits)) && (archDescriptionList$targetArch=="singleCPU")) {
        while (length(activeExpUnits)) {
            currentLevel <- activeExpUnits[1]
            activeExpUnits <- activeExpUnits[-1]
            .clusterExpUnitWithScamp(
                expUnit=currentLevel,
                resList=resList,
                selectedChannels=selectedChannels,
                numScampIter=numScampIter,
                debugFlag=debugFlag,
                threadNum=threadNum,
                seedValue=seedValue,
                densitySubSampleThreshold=densitySubSampleThreshold,
                densitySubSampleSize=densitySubSampleSize,
                densitySubSampleIterations=densitySubSampleIterations,
                projectPath=projectPath
            )
        }
    }
    else if ((length(activeExpUnits)) && (archDescriptionList$targetArch=="slurmCluster")) {
        if (!dir.exists(file.path(normalizePath(projectPath),
                                  "faustData",
                                  "slurmScampData"))) {
            dir.create(file.path(normalizePath(projectPath),
                                 "faustData",
                                 "slurmScampData"))
        }
        stillRunningSlurm <- TRUE
        startSlurmTime <- proc.time()
        maxNodeNum <- archDescriptionList$maxNodeNum
        maxTime <- archDescriptionList$maxTime
        sbatchFlags <- archDescriptionList$sbatchFlags
        currentJobs <- 0
        jobNum <- 0
        slurmLevels <- c()
        while (stillRunningSlurm) {
            if ((currentJobs < maxNodeNum) && (length(activeExpUnits))) {
                jobNum <- jobNum + 1
                currentLevel <- activeExpUnits[1]
                activeExpUnits <- activeExpUnits[-1]
                currentJobs <- (currentJobs + 1)
                slurmLevels <- append(slurmLevels,currentLevel)
                if (!dir.exists(file.path(
                         normalizePath(projectPath),
                         "faustData",
                         "slurmScampData",
                         currentLevel
                     )))
                {
                    dir.create(file.path(normalizePath(projectPath),
                                         "faustData",
                                         "slurmScampData",
                                         currentLevel))
                }
                .prepareSlurmScampJob(
                    expUnit = currentLevel,
                    startingCellPop = startingCellPop,
                    selectedChannels = selectedChannels,
                    numScampIter = numScampIter,
                    minClusterSize = 25,
                    threadNum = archDescriptionList$nodeThreadNum,
                    seedValue = seedValue,
                    projectPath = projectPath,
                    jobNumber = jobNum,
                    partitionID = archDescriptionList$partitionID,
                    jobTime = archDescriptionList$jobTime,
                    jobPrefix = archDescriptionList$jobPrefix,
                    densitySubSampleThreshold = densitySubSampleThreshold,
                    densitySubSampleSize = densitySubSampleSize,
                    densitySubSampleIterations = densitySubSampleIterations
                )
                print(paste0("Slurm SCAMP clustering starting for ", currentLevel))
                launchJob <- system2("sbatch",
                                     args=paste0(sbatchFlags,
                                                 paste0(" '",
                                                        file.path(
                                                            normalizePath(projectPath),
                                                            "faustData",
                                                            "slurmScampData",
                                                            currentLevel,
                                                            "slurmScampJob.sh"
                                                        ),
                                                        "'")
                                                 ),
                                     stdout=TRUE)
            }
            else {
                Sys.sleep(10) #in seconds
                currentSlurmTime <- (proc.time() - startSlurmTime)
                if (as.numeric(currentSlurmTime[3]) > maxTime) {
                    print("Slurm SCAMP clustering exceeded max time.")
                    print(paste0("Check logs in ",file.path("faustData","slurmScampData")))
                    stop("Killing FAUST")
                }
                activeSlurmLevels <- c()
                for (sLevel in slurmLevels) {
                    if ((file.exists(file.path(normalizePath(projectPath),
                                               "faustData",
                                               "slurmScampData",
                                               sLevel,
                                               "slurmScampComplete.rds"))) &&
                        (readRDS(file.path(normalizePath(projectPath),
                                           "faustData",
                                           "slurmScampData",
                                           sLevel,
                                           "slurmScampComplete.rds"))))
                    {
                        print(paste0("Slurm SCAMP clustering complete for ", sLevel))
                        currentJobs <- (currentJobs - 1)
                    }
                    else {
                        activeSlurmLevels <- append(activeSlurmLevels,sLevel)
                    }
                }
                slurmLevels <- activeSlurmLevels
                if ((length(activeExpUnits)==0) && (currentJobs==0)) {
                    stillRunningSlurm <- FALSE
                }
            }
        }
    }
    else if (length(activeExpUnits)) {
        print("Unsupported targetArch requested in archDescriptionList.")
        stop("Killing FAUST.")
    }
    #collect all labels from the scamp clusterings
    if (debugFlag) print("Accumulating cluster labels.")
    clusterNames <- c()
    for (experimentalUnit in uniqueExpUnits) {
        if (!file.exists(file.path(normalizePath(projectPath),
                                   "faustData",
                                   "expUnitData",
                                   experimentalUnit,
                                   "scampClusterLabels.rds"))) {
            print(paste0("Labels not detected in experimental unit ",experimentalUnit))
            print("This is a bug -- all experimental units should have labels.")
            stop("Killing FAUST. Check logs to determine which experimental unit is unlabeled.")
        }
        else {
            expUnitLabels <- readRDS(file.path(normalizePath(projectPath),
                                             "faustData",
                                             "expUnitData",
                                             experimentalUnit,
                                             "scampClusterLabels.rds"))
            clusterNames <- append(clusterNames,expUnitLabels)
        }
    }
    nameSummary <- table(clusterNames)
    saveRDS(nameSummary,
            file.path(normalizePath(projectPath),
                      "faustData",
                      "metaData",
                      "scampNameSummary.rds"))
    nameSummaryPlotDF <- data.frame(x=seq(max(nameSummary)),
                                    y=sapply(seq(max(nameSummary)),function(x){
                                        length(which(nameSummary >= x))}))
    
    elbowLoc <- .computeElbow(nameSummaryPlotDF)
    if (nameOccurenceNum > 0) {
        #the user has set this value
        clusterNames <- names(nameSummary[which(nameSummary >= nameOccurenceNum)])
        saveRDS(
            nameOccurenceNum,
            file.path(normalizePath(projectPath),
                      "faustData",
                      "metaData",
                      "phenotypeElbowValue.rds")
        )
    }
    else {
        #use the automatic value.
        clusterNames <- names(nameSummary[which(nameSummary >= elbowLoc)])
        saveRDS(
            elbowLoc,
            file.path(normalizePath(projectPath),
                      "faustData",
                      "metaData",
                      "phenotypeElbowValue.rds")
        )
    }
    saveRDS(clusterNames,
            file.path(normalizePath(projectPath),
                      "faustData",
                      "metaData",
                      "scampClusterNames.rds"))

    allClusterNames <- names(nameSummary[which(nameSummary >= 1)])
    saveRDS(allClusterNames,
            file.path(normalizePath(projectPath),
                      "faustData",
                      "metaData",
                      "allScampClusterNames.rds"))
    
    if (debugFlag) print("Cluster labels collected and saved.")
    return()
}



.prepareSlurmScampJob <- function(expUnit,
                                  startingCellPop,
                                  selectedChannels,
                                  numScampIter,
                                  minClusterSize,
                                  debugFlag,
                                  threadNum,
                                  seedValue,
                                  projectPath,
                                  jobNumber,
                                  partitionID,
                                  jobTime,
                                  jobPrefix,
                                  densitySubSampleThreshold,
                                  densitySubSampleSize,
                                  densitySubSampleIterations
                                  )
{
    .programTemplate <-'library(scamp)
library(data.table)
expUnitExprs <- readRDS(file.path(normalizePath({{projectPath}}),"faustData","expUnitData",{{expUnit}},"expUnitExprs.rds"))
expUnitRes <- readRDS(file.path(normalizePath({{projectPath}}),"faustData","expUnitData",{{expUnit}},"expUnitRes.rds"))
resList <- readRDS(file.path(normalizePath({{projectPath}}),"faustData","gateData",paste0({{startingCellPop}},"_resList.rds")))
expUnitExprs <- expUnitExprs[,{{selectedChannels}}, drop = FALSE]
expUnitRes <- expUnitRes[,{{selectedChannels}}, drop = FALSE]
resFlag <- FALSE
for (colNum in seq(ncol(expUnitRes))) {
    if (length(which(expUnitRes[,colNum] > 0))) {
        resFlag <- TRUE
        break
    }
}
#get the scamp annotation forest
scampAF <- rep(list(NA),ncol(expUnitExprs))
names(scampAF) <- colnames(expUnitExprs)
for (colName in colnames(expUnitExprs)) {
    scampAF[[colName]] <- resList[[colName]][[{{expUnit}}]]
}
scampClustering <- scamp::scamp(
    dataSet = expUnitExprs,
    numberIterations = {{numScampIter}},
    minimumClusterSize = {{minClusterSize}},
    numberOfThreads = {{threadNum}},
    anyValueRestricted = resFlag,
    resValMatrix = expUnitRes,
    useAnnForest = TRUE,
    annForestVals = scampAF,
    randomSeed={{seedValue}},
    getDebugInfo = FALSE,
    subSampleThreshold = {{densitySubSampleThreshold}},
    subSampleSize = {{densitySubSampleSize}},
    subSampleIterations = {{densitySubSampleIterations}}
)
runClustering <- scampClustering[[1]]
maxClustering <- scampClustering[[2]]
outClustering <- rep("Uncertain",length(maxClustering))
agreeIndex <- which(runClustering==maxClustering)
if (length(agreeIndex) == 0) {
   print(paste0("SCAMP reached no concensus for: ",expUnit))
   print("This indicates high-amount of uncertainty in cluster labels.")
   print("Increase the numScampIter parameters until concensus reached")
   stop("killing job.")
}
outClustering[agreeIndex] <- maxClustering[agreeIndex]
clusterNames <- setdiff(sort(unique(names(table(outClustering)))),"Uncertain")
saveRDS(clusterNames,file.path(normalizePath({{projectPath}}),"faustData","expUnitData",{{expUnit}},"scampClusterLabels.rds"))
#unwind the experimental unit to each sample
expUnitToSampleLookup <- readRDS(file.path(normalizePath({{projectPath}}),"faustData","expUnitData",{{expUnit}},"expUnitToSampleLookup.rds"))
for (sampleName in names(table(expUnitToSampleLookup))) {
    sampleLookup <- which(expUnitToSampleLookup == sampleName)
    if (length(sampleLookup)) {
        sampleClustering <- outClustering[sampleLookup]
        data.table::fwrite(
            list(sampleClustering),
            file = file.path(normalizePath({{projectPath}}),"faustData","sampleData",sampleName,"scampAnnotation.csv"),
            sep = "`",
            append = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)
    }
}
scampExpUnitDone <- TRUE
saveRDS(scampExpUnitDone,file.path(normalizePath({{projectPath}}),"faustData","expUnitData",{{expUnit}},"scampExpUnitComplete.rds"))
slurmScampDone <- TRUE
saveRDS(slurmScampDone,file.path(normalizePath({{projectPath}}),"faustData","slurmScampData",{{expUnit}},"slurmScampComplete.rds"))
'
    programData <- list(
        expUnit=paste0("'",expUnit,"'"),
        startingCellPop=paste0("'",startingCellPop,"'"),
        selectedChannels=paste0("c('",paste0(selectedChannels,collapse="','"),"')"),
        numScampIter=numScampIter,
        minClusterSize=minClusterSize,
        threadNum=threadNum,
        seedValue=seedValue,
        projectPath=paste0("'",projectPath,"'"),
        densitySubSampleThreshold=densitySubSampleThreshold,
        densitySubSampleSize=densitySubSampleSize,
        densitySubSampleIterations=densitySubSampleIterations
    )
    renderedProgram <- whisker.render(.programTemplate, programData)
    write(
        renderedProgram,
        file=file.path(normalizePath(projectPath),
                       "faustData",
                       "slurmScampData",
                       expUnit,
                       "slurmScampJob.R")
    )
    .controlTemplate <-'#!/bin/bash
#SBATCH --partition={{partitionID}}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={{threadNum}}
#SBATCH --time={{jobTime}}
#SBATCH -o {{logPath}}
#SBATCH -J {{jobPrefix}}sj{{jobNumber}}
#SBATCH --threads-per-core=1

echo "Start of program at `date`"
Rscript --no-save --no-restore {{jobPath}}
echo "End of program at `date`"'
    controlData <- list(
        threadNum = threadNum,
        jobNumber = jobNumber,
        partitionID = partitionID,
        jobTime = jobTime,
        jobPrefix = jobPrefix,
        jobPath = paste0("'",
                         file.path(normalizePath(projectPath),
                                   "faustData",
                                   "slurmScampData",
                                   expUnit,
                                   "slurmScampJob.R"),
                         "'"),
        logPath = paste0("'",
                         file.path(normalizePath(projectPath),
                                   "faustData",
                                   "slurmScampData",
                                   expUnit,
                                   "sjLog"),
                         "'")
    )
    renderedScript <- whisker.render(.controlTemplate, controlData)
    write(
        renderedScript,
        file=file.path(normalizePath(projectPath),
                       "faustData",
                       "slurmScampData",
                       expUnit,
                       "slurmScampJob.sh")
    )
    return()
}
