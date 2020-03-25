.growForestForExpUnit <- function(expUnit,
                                 rootPop,
                                 activeChannels,
                                 analysisMap,
                                 numIter,
                                 debugFlag,
                                 threadNum,
                                 seedValue,
                                 projectPath,
                                 densitySubSampleThreshold,
                                 densitySubSampleSize,
                                 densitySubSampleIterations)
{
    #function to
    if (debugFlag) print(paste0("Growing annotation forest for: ",expUnit))
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
    expUnitExprs <- expUnitExprs[,activeChannels,drop=FALSE]
    expUnitRes <- expUnitRes[,activeChannels,drop=FALSE]
    resFlag <- FALSE
    for (colNum in seq(ncol(expUnitRes))) {
        if (length(which(expUnitRes[,colNum] > 0))) {
            resFlag <- TRUE
            break
        }
    }
    annF <- growAnnotationForest(
        dataSet=expUnitExprs,
        numberIterations=numIter,
        pValueThreshold=0.25,
        minimumClusterSize=25,
        randomCandidateSearch=FALSE,
        maximumSearchDepth=2,
        numberOfThreads=threadNum,
        maximumGatingNum=1e10,
        anyValueRestricted=resFlag,
        resValMatrix=expUnitRes,
        cutPointUpperBound=2,
        getDebugInfo=FALSE,
        randomSeed=seedValue,
        subSamplingThreshold=densitySubSampleThreshold,
        subSampleSize=densitySubSampleSize,
        subSampleIter=densitySubSampleIterations,
        recordCounts=FALSE,
        recordIndices=FALSE
    )
    saveRDS(annF,
            file.path(normalizePath(projectPath),
                      "faustData",
                      "expUnitData",
                      expUnit,
                      paste0(rootPop,"_annF.rds")))
    ePop <- apply(expUnitRes,2,function(x){length(which(x==0))})
    names(ePop) <- colnames(expUnitExprs)
    af <- annF[["gateData"]]
    pAnnF <- .parseAnnotationForest(af,ePop)
    saveRDS(pAnnF,
            file.path(normalizePath(projectPath),
                      "faustData",
                      "expUnitData",
                      expUnit,
                      paste0(rootPop,"_pAnnF.rds")))
    expUnitDone <- TRUE
    saveRDS(expUnitDone,
            file.path(normalizePath(projectPath),
                      "faustData",
                      "expUnitData",
                      expUnit,
                      "expUnitComplete.rds"))
    if (debugFlag) print(paste0("Annotation forest complete for: ",expUnit))
    return()
}



.growAnnForest <- function(activeChannels,
                           debugFlag,
                           threadNum,
                           seedValue,
                           projectPath,
                           densitySubSampleThreshold,
                           densitySubSampleSize,
                           densitySubSampleIterations,
                           archDescriptionList)
{
    #removing numForestIter from interface for simplicity,
    #since it is a rarely modified parameter that makes
    #it harder to apply faust for a new user.
    #if we end up wanting to add numForestIter back to the
    #interface, exposed the numIter parameter here.
    numIter <- 1
    rootPop <- readRDS(file.path(normalizePath(projectPath),
                                 "faustData",
                                 "metaData",
                                 "sanitizedCellPopStr.rds"))

    analysisMap <- readRDS(file.path(normalizePath(projectPath),
                                     "faustData",
                                     "metaData",
                                     "analysisMap.rds"))

    uniqueExpUnits <- sort(unique(analysisMap[,"experimentalUnit"]))
    activeExpUnits <- c()
    #accumulate vector of experimental units without annotation forests.
    for (experimentalUnit in uniqueExpUnits) {
        if (!file.exists(file.path(normalizePath(projectPath),
                                   "faustData",
                                   "expUnitData",
                                   experimentalUnit,
                                   "expUnitComplete.rds"))) {
            activeExpUnits <- append(activeExpUnits,experimentalUnit)
        }
    }
    #grow forests for experimental units that lack them
    if ((length(activeExpUnits)) && (archDescriptionList$targetArch=="singleCPU")) {
        while (length(activeExpUnits)) {
            currentLevel <- activeExpUnits[1]
            activeExpUnits <- activeExpUnits[-1]
            .growForestForExpUnit(
                expUnit=currentLevel,
                rootPop=rootPop,
                activeChannels=activeChannels,
                analysisMap=analysisMap,
                numIter=numIter,
                debugFlag=debugFlag,
                threadNum=threadNum,
                seedValue=seedValue,
                projectPath=projectPath,
                densitySubSampleThreshold=densitySubSampleThreshold,
                densitySubSampleSize=densitySubSampleSize,
                densitySubSampleIterations=densitySubSampleIterations
            )
        }
    }
    else if ((length(activeExpUnits)) && (archDescriptionList$targetArch=="slurmCluster")) {
        if (!dir.exists(file.path(normalizePath(projectPath),
                                  "faustData",
                                  "slurmData"))) {
            dir.create(file.path(normalizePath(projectPath),
                                 "faustData",
                                 "slurmData"))
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
                if (!dir.exists(file.path(normalizePath(projectPath),
                                          "faustData",
                                          "slurmData",
                                          currentLevel))) {
                    dir.create(file.path(normalizePath(projectPath),
                                         "faustData",
                                         "slurmData",
                                         currentLevel))
                }
                .prepareSlurmJob(
                    expUnit = currentLevel,
                    rootPop = rootPop,
                    activeChannels = activeChannels,
                    analysisMap = analysisMap,
                    numIter = numIter,
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
                print(paste0("Slurm annotation forest starting for ", currentLevel))
                launchJob <- system2("sbatch",
                                     args=paste0(sbatchFlags,
                                                 paste0(" '",
                                                        file.path(
                                                            normalizePath(projectPath),
                                                            "faustData",
                                                            "slurmData",
                                                            currentLevel,
                                                            "slurmJob.sh"
                                                        ),
                                                        "'")
                                                 ),
                                     stdout=TRUE)
            }
            else {
                Sys.sleep(10) #in seconds
                currentSlurmTime <- (proc.time() - startSlurmTime)
                if (as.numeric(currentSlurmTime[3]) > maxTime) {
                    print("Slurm annotation forest exceeded max time.")
                    print(paste0("Check logs in ",
                                 file.path(normalizePath(projectPath),
                                           "faustData",
                                           "slurmData")))
                    stop("Killing FAUST")
                }
                activeSlurmLevels <- c()
                for (sLevel in slurmLevels) {
                    if (
                    (file.exists(file.path(normalizePath(projectPath),
                                           "faustData",
                                           "slurmData",
                                           sLevel,
                                           "slurmComplete.rds"))) &&
                    (readRDS(file.path(normalizePath(projectPath),
                                       "faustData",
                                       "slurmData",
                                       sLevel,
                                       "slurmComplete.rds")))
                    )
                    {
                        print(paste0("Slurm annotation forest complete for ", sLevel))
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
    return()
}

.prepareSlurmJob <- function(expUnit,
                             rootPop,
                             activeChannels,
                             analysisMap,
                             numIter,
                             threadNum,
                             seedValue,
                             projectPath,
                             jobNumber,
                             partitionID,
                             jobTime,
                             jobPrefix,
                             densitySubSampleThreshold,
                             densitySubSampleSize,
                             densitySubSampleIterations)
{
    .programTemplate <-'library(faust)
expUnitExprs <- readRDS(file.path(normalizePath({{projectPath}}),"faustData","expUnitData",{{expUnit}},"expUnitExprs.rds"))
expUnitRes <- readRDS(file.path(normalizePath({{projectPath}}),"faustData","expUnitData",{{expUnit}},"expUnitRes.rds"))
expUnitExprs <- expUnitExprs[,{{activeChannels}},drop=FALSE]
expUnitRes <- expUnitRes[,{{activeChannels}},drop=FALSE]
resFlag <- FALSE
for (colNum in seq(ncol(expUnitRes))) {
if (length(which(expUnitRes[,colNum] > 0))) {
resFlag <- TRUE
break
}
}
annF <- faust:::growAnnotationForest(
dataSet=expUnitExprs,
numberIterations={{numIter}},
pValueThreshold=0.25,
minimumClusterSize=25,
randomCandidateSearch=FALSE,
maximumSearchDepth=2,
numberOfThreads={{threadNum}},
maximumGatingNum=1e10,
anyValueRestricted=resFlag,
resValMatrix=expUnitRes,
cutPointUpperBound=2,
getDebugInfo=FALSE,
randomSeed={{seedValue}},
subSamplingThreshold={{densitySubSampleThreshold}},
subSampleSize={{densitySubSampleSize}},
subSampleIter={{densitySubSampleIterations}},
recordCounts=FALSE,
recordIndices=FALSE
)
saveRDS(annF,file.path(normalizePath({{projectPath}}),"faustData","expUnitData",{{expUnit}},paste0({{rootPop}},"_annF.rds")))
ePop <- apply(expUnitRes,2,function(x){length(which(x==0))})
names(ePop) <- colnames(expUnitExprs)
af <- annF[["gateData"]]
pAnnF <- faust:::.parseAnnotationForest(af,ePop)
saveRDS(pAnnF,file.path(normalizePath({{projectPath}}),"faustData","expUnitData",{{expUnit}},paste0({{rootPop}},"_pAnnF.rds")))
expUnitDone <- TRUE
saveRDS(expUnitDone,file.path(normalizePath({{projectPath}}),"faustData","expUnitData",{{expUnit}},"expUnitComplete.rds"))
slurmDone <- TRUE
saveRDS(slurmDone,file.path(normalizePath({{projectPath}}),"faustData","slurmData",{{expUnit}},"slurmComplete.rds"))
'
    programData <- list(
        expUnit=paste0("'",expUnit,"'"),
        rootPop=paste0("'",rootPop,"'"),
        activeChannels=paste0("c('",paste0(activeChannels,collapse="','"),"')"),
        numIter=numIter,
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
                       "slurmData",
                       expUnit,
                       "slurmJob.R")
    )
    .controlTemplate <-'#!/bin/bash
#SBATCH --partition={{partitionID}}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={{threadNum}}
#SBATCH --time={{jobTime}}
#SBATCH -o {{logPath}}
#SBATCH -J {{jobPrefix}}fj{{jobNumber}}
#SBATCH --threads-per-core=1

echo "Start of program at `date`"
Rscript --no-save --no-restore {{jobPath}}
echo "End of program at `date`"'
    controlData <- list(
        threadNum = threadNum,
        jobNumber = jobNumber,
        partitionID = partitionID,
        jobTime=jobTime,
        jobPrefix=jobPrefix,
        jobPath = paste0("'",
                         file.path(normalizePath(projectPath),
                                   "faustData",
                                   "slurmData",
                                   expUnit,
                                   "slurmJob.R"),
                         "'"),
        logPath = paste0("'",
                         file.path(normalizePath(projectPath),
                                   "faustData",
                                   "slurmData",
                                   expUnit,
                                   "fjLog"),
                         "'")
    )
    renderedScript <- whisker.render(.controlTemplate, controlData)
    write(
        renderedScript,
        file=file.path(normalizePath(projectPath),
                       "faustData",
                       "slurmData",
                       expUnit,
                       "slurmJob.sh")
    )
    return()
}

