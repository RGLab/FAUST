library(FlowSOM)
library(lme4)
library(tidyr)
library(dplyr)
library(faust)
library(MCMCpack)
library(mvtnorm)
library(MASS)
library(openCyto)
library(flowCore)
library(ggplot2)
library(optparse)
library(gmp)

option_list = list(make_option(c("-a", "--numOfClusters"),
                               action="store",
                               default=NA,
                               type="character",
                               help="Subject to process"),
                   make_option(c("-b", "--noiseNumber"),
                               action="store",
                               default=NA,
                               type="character",
                               help="Subject to process"),
                   make_option(c("-c", "--batchNumber"),
                               action="store",
                               default=NA,
                               type="character",
                               help="Subject to process"),
                   make_option(c("-d", "--transNumber"),
                               action="store",
                               default=NA,
                               type="character",
                               help="Subject to process"))

opt = parse_args(OptionParser(option_list=option_list))

#get simulation choices from the command line.
clusterChoice <- as.numeric(opt$numOfClusters)
noiseChoice <- as.numeric(opt$noiseNumber)
batchChoice <- as.numeric(opt$batchNumber)
transformationChoice <- as.numeric(opt$transNumber)

#uncomment the following lines if you want to interactively set simulation parameters.
#clusterChoice <- 1
#batchChoice <- 1
#noiseChoice <- 1
#transformationChoice <- 1

#possible settings.
possibleClusters <- c(10,25,50)
batchEffects <- c(0,(1/3))
noiseDimensions <- c(0,5)
possibleOPGridNums <- c(5,7,10)

#the dimension of spec mat determines the true annotation of clusters in the simulated experiment.
#a mean vector for the multivariate normal is sampled from each column
#if the 0 is selected, the true annotation is "low" for that variable.
#if the non-0 is selected, the true annotation is "high" for that variable.
specMat <- matrix(0,nrow=2,ncol=10)
specMat[2,] <- c(rep(8,2),rep(7,2),rep(6,2),rep(5,2),rep(4,2))
if (transformationChoice == 1) {
    dataTransformations <- list(function(x){return(x)},
                                function(x){return(x)},
                                function(x){return(x)})
}
if (transformationChoice == 2) {
    dataTransformations <- list(function(x){return((x^2))},
                                function(x){return((x^2))},
                                function(x){return((x^2))})
}
if (transformationChoice == 3) {
    #adjust specification for cytof-like transformation
    specMat[2,] <- c(rep(8,3),rep(7,4),rep(6,3))
    dataTransformations <- list(function(x){return(gamma(1+abs(x/4)))},
                                function(x){return(gamma(1+abs(x/4)))},
                                function(x){return(gamma(1+abs(x/4)))})
}
print("specMat")
print(specMat)
print("dataTransformations")
print(dataTransformations)
#select simulation settings
numberOfClusters <- possibleClusters[clusterChoice]
opGridNum <- possibleOPGridNums[clusterChoice]
#
#the "spiked in" population has increased abundance in the spiked setting
#in the unspiked setting, the population is the smallest (out of number of clusters)
#
if (numberOfClusters==10) {
    targetSpikeRank <- 7
}
if (numberOfClusters==25) {
    targetSpikeRank <- 18
}
if (numberOfClusters==50) {
    targetSpikeRank <- 35
}
print(paste0("Targeting spike: ",targetSpikeRank))
batchEffectSize <- batchEffects[batchChoice]
#scale batch effect size to transfromation type -- the batch is incremented 9 times over the simulation.
if ((batchEffectSize > 0) && (transformationChoice == 2)) {
    batchEffectSize <- (1/3)
}
if ((batchEffectSize > 0) && (transformationChoice == 3)) {
    batchEffectSize <- 0.04923483
}
nuisanceVariables <- noiseDimensions[noiseChoice]

#load functions for simulation
source("./functionsForSimulation.R")

#report settings to log.
print(paste0("Selected data transformations: ", dataTransformations))
print(paste0("Selected number of clusters: ", numberOfClusters))
print(paste0("Selected batch effect: ", batchEffectSize))
print(paste0("Selected nuisance: ", nuisanceVariables))
print(paste0("Selected over partitioned grid: ", opGridNum))
jobIdStr <- paste0(
    "c",clusterChoice,
    "b",batchChoice,
    "n",noiseChoice,
    "t",transformationChoice
)
jobPath <- paste0("./simResults/",jobIdStr)
fileAbortPath <- paste0(jobPath,"/stopJob.txt")
resultsPath <- paste0(jobPath,"/simOutput.rds")
tmpResultsPath <- paste0(jobPath,"/intermediateSimOutput.rds")
faustDataPath <- paste0(jobPath,"/faustData")
print(paste0("jobPath: ",jobPath))
print(paste0("fileAbortPath: ",fileAbortPath))
print(paste0("resultsPath: ",resultsPath))
print(paste0("faustDataPath: ",faustDataPath))

#create directories for reporting
if (!dir.exists(jobPath)) {
    dir.create(jobPath)
}

#
#set seed and start the simulation.
#
rsValue <- 777
print(paste0("Random seed to start simulation: ", rsValue))
set.seed(rsValue)
allSimulationResults <- list()
#
#modify simIter to control the number of total iterations run for the sim settings.
#
simIter <- 50
while(simIter) {
print(paste0("beginning sim iter: ",simIter))
#
#begin by sampling a probability vector across the number of clusters.
#note that the final cluster is always the spiked in element (hence re-arrangement)
#in the base scenario, set the spiked in element to have minimum mass.
#
minVal <- 0
while (minVal < 1e-3) {
    clusterProbVec <- rdirichlet(1, rep(1,numberOfClusters))
    minVal <- min(clusterProbVec)
}
mvLookup <- which(clusterProbVec == minVal)[1]
lastMass <- clusterProbVec[length(clusterProbVec)]
if (mvLookup != length(clusterProbVec)) {
    clusterProbVec[mvLookup] <- lastMass
    clusterProbVec[length(clusterProbVec)] <- minVal
}
if (abs(sum(clusterProbVec)-1) > 1e-9) {
    stop("Initial probability error.")
}
if (min(clusterProbVec) != clusterProbVec[length(clusterProbVec)]) {
    stop("min elt not final elt")
}
                 
#
#Simulate an experiment.
#
simmedExp <- simulateExperiment(
    meanVectorBoundsMatrix=specMat,
    numSamples=100,
    transformationList=dataTransformations,
    batchEffectShift=batchEffectSize,
    noiseDimension=nuisanceVariables,
    probVecIn=clusterProbVec,
    minSampleSize=5000,
    tncp=10000,
    targetRank=targetSpikeRank
)
#
#uncomment to see the distribution of spiked in population stratified by responder status
#
#chkDF <- data.frame(pop=simmedExp[["truthMat"]][,simmedExp[["spikedPop"]]],
#                    resp=as.factor(unlist(simmedExp[["spikedInPopList"]])))
#ggplot(chkDF,aes(x=resp,y=pop))+geom_boxplot()+theme_bw()

#
#discover and annotate the data with FAUST
#
activeChannelsIn <- paste0("V",seq((10 + nuisanceVariables)))
print("activeChannelsIn")
print(activeChannelsIn)
#channelBoundsIn <- matrix(-Inf,nrow=2,ncol=length(activeChannelsIn))
#colnames(channelBoundsIn) <- activeChannelsIn
#rownames(channelBoundsIn) <- c("Low","High")
#channelBoundsIn["High",] <- Inf
fs <- as(simmedExp[["flowFrameList"]], "flowSet")
gs <- GatingSet(fs)
faust:::faust(
    gatingSet=gs,
    activeChannels=activeChannelsIn,
    #channelBounds=channelBoundsIn,
    startingCellPop="root",
    projectPath=jobPath,
    experimentalUnit="name",
    depthScoreThreshold = 0.01,
    selectionQuantile = 0.5,
    drawAnnotationHistograms = FALSE,
    nameOccurenceNum = max(1,floor(length(gs)/4)),
    debugFlag = TRUE,
    threadNum = 10,
    seedValue = 12345,
    annotationsApproved = TRUE,
    plottingDevice = "pdf"
)
#
#Concatenate and cluster the data with flowSOM
#
totalL <- sum(unlist(lapply(simmedExp[["flowFrameList"]],nrow)))
trueLabelVec <- sampleLookups  <- rep("NoSample",totalL)
exprsMat <- matrix(0,nrow=totalL,ncol=length(activeChannelsIn))
colnames(exprsMat) <- activeChannelsIn
startL <- 1
for (sampleName in sampleNames(gs)) {
    print(sampleName)
    newData <- exprs(getData(gs[[sampleName]],"root"))
    newData <- newData[,activeChannelsIn]
    endL <- (startL + nrow(newData) - 1)
    exprsMat[startL:endL,] <- newData
    sampleLookups[startL:endL] <- sampleName
    trueLabelVec[startL:endL] <- simmedExp[["labelsList"]][[which(names(simmedExp[["labelsList"]])==sampleName)]]
    startL <- (endL + 1)
}
ffFlowSOM <- flowCore::flowFrame(exprsMat)
fSOM <- FlowSOM::ReadInput(ffFlowSOM, transform = FALSE, scale = TRUE)
#
#
#First, we set the grid to 1 by number of true clusters (mixture components).
#The choice of the 1 by number of true clusters comes from
#Handbook of Cluster Analysis, Edited by Christian Hennig, Marina Meila
#Fionn Murtagh, and Roberto Rocci, c 2016 Taylor & Francis Group LLC 
#Version Date: 20151012, ISBN-13: 978-1-4665-5188-6
#SOM discussed pages 421-423. Clustering choice page 421.
#
#
fSOM <- FlowSOM::BuildSOM(fSOM,
                          colsToUse = seq(length(colnames(exprsMat))),
                          xdim = 1,
                          ydim = length(clusterProbVec))
fSOM <- FlowSOM::BuildMST(fSOM)
clusterLabels <- fSOM$map$mapping[, 1]
#
#record those clusters containing the maximum spiked in population and the highest % spiked in 
#
modOCSummary <- table(clusterLabels)
modOCMSummary <- table(clusterLabels[which(trueLabelVec==simmedExp[["spikedPop"]])])
modMaxCountOracleNames <- names(which(modOCMSummary==max(modOCMSummary)))
modOracleProps <- c()
for (clusterName in names(modOCMSummary)) {
    modOracleProps <- append(modOracleProps,(modOCMSummary[clusterName]/modOCSummary[clusterName]))
}
modMaxPropOracleNames <- names(which(modOracleProps==max(modOracleProps)))
#
#derive the count matrix by sample from the flowSOM clustering.
#
sampleNames <- names(table(sampleLookups))
flowsomCountMatrix <- matrix(0,nrow=length(sampleNames),ncol=length(clusterProbVec))
rownames(flowsomCountMatrix) <- sampleNames
colnames(flowsomCountMatrix) <- seq(length(clusterProbVec))
for (sampleName in sampleNames) {
    rowLookup <- which(rownames(flowsomCountMatrix)==sampleName)
    dataLookup <- which(sampleLookups==sampleName)
    sampleCluster <- clusterLabels[dataLookup]
    for (clusterNum in seq(length(clusterProbVec))) {
        clusterCount <- length(which(sampleCluster==clusterNum))
        if (clusterCount > 0) {
            colLookup <- which(colnames(flowsomCountMatrix)==as.character(clusterNum))
            flowsomCountMatrix[rowLookup,colLookup] <- clusterCount
        }
    }
}

#
#
#We next set a grid that is deliberately over-partitioned, and cluster the data
#with FlowSOM again.
#
#
opfSOM <- FlowSOM::BuildSOM(fSOM,
                            colsToUse = seq(length(colnames(exprsMat))),
                            xdim = opGridNum,
                            ydim = opGridNum)
opfSOM <- FlowSOM::BuildMST(opfSOM)
opClusterLabels <- opfSOM$map$mapping[, 1]
#
#record those clusters containing the maximum spiked in population and the highest % spiked in 
#
modOPSummary <- table(opClusterLabels)
modOPMSummary <- table(opClusterLabels[which(trueLabelVec==simmedExp[["spikedPop"]])])
modMaxCountOverpartNames <- names(which(modOPMSummary==max(modOPMSummary)))
modOverpartProps <- c()
for (clusterName in names(modOPMSummary)) {
    modOverpartProps <- append(modOverpartProps,(modOPMSummary[clusterName]/modOPSummary[clusterName]))
}
modMaxPropOverpartNames <- names(which(modOverpartProps==max(modOverpartProps)))

#
#Next, derive the count matrix by sample from the overpartitioned flowSOM clustering.
#
overpartitionedCountMatrix <- matrix(0,nrow=length(sampleNames),ncol=(opGridNum * opGridNum))
rownames(overpartitionedCountMatrix) <- sampleNames
colnames(overpartitionedCountMatrix) <- seq((opGridNum*opGridNum))
for (sampleName in sampleNames) {
    rowLookup <- which(rownames(overpartitionedCountMatrix)==sampleName)
    dataLookup <- which(sampleLookups==sampleName)
    sampleCluster <- opClusterLabels[dataLookup]
    for (clusterNum in seq((opGridNum*opGridNum))) {
        clusterCount <- length(which(sampleCluster==clusterNum))
        if (clusterCount > 0) {
            colLookup <- which(colnames(overpartitionedCountMatrix)==as.character(clusterNum))
            overpartitionedCountMatrix[rowLookup,colLookup] <- clusterCount
        }
    }
}

#
#Check to make sure the counts match by sample for the two FlowSOM clusterings, but columns differ.
#
flowsomCountsMatch <- min(apply(overpartitionedCountMatrix,1,sum)==apply(flowsomCountMatrix,1,sum))
if (!flowsomCountsMatch) {
    stop("Error in deriving overparititoned grid.")
}
if (ncol(flowsomCountMatrix)==ncol(overpartitionedCountMatrix)) {
    stop("Error in deriving overparititoned grid.")
}

#
#Read in the Faust counts.
#
faustMat <- readRDS(paste0(faustDataPath,"/faustCountMatrix.rds"))
medianFaustPctUnclassified <- median(faustMat[,"0_0_0_0_0"]/apply(faustMat,1,sum))

#
#Reorder the faust results and the true annotations to follow depth-score order (for comparability)
#
fo <- getStrOrder(colnames(faustMat)[1])
faustVarSelectionCount <- length(which(fo <= 10))
newLabels <- as.character(sapply(colnames(simmedExp[["truthMat"]]),function(x){truth2faust(x,fo)}))
modTruthMat <- simmedExp[["truthMat"]]
colnames(modTruthMat) <- newLabels
pcTrue <- apply(modTruthMat,1,sum)
trueProps <- apply(modTruthMat,2,function(x){x/pcTrue})
allSubjects <- rownames(modTruthMat)
spikedInPop <- truth2faust(simmedExp[["spikedPop"]],fo)

#
#Record if FAUST finds the true populations, how it correlates with them.
#
matchingNamesLookup <- which(colnames(faustMat) %in% colnames(modTruthMat))
spikedCorrelation <- medFaustCorrelation <- spikedInFound <- truePopsFound <- 0
if (length(matchingNamesLookup)) {
    matchingNames <- colnames(faustMat)[matchingNamesLookup]
    truePopsFound <- length(matchingNames)
    corVec <- c()
    for (cName in matchingNames) {
        corVec <- append(corVec,cor(modTruthMat[,cName],faustMat[,cName]))
    }
    medFaustCorrelation <- median(corVec)
    spikeLookup <- which(matchingNames==spikedInPop)
    if (length(spikeLookup)) {
        spikedInFound <- 1
        spikedCorrelation <- corVec[spikeLookup]
    }
}
faustDF <- as.data.frame(faustMat)
cellPops <- colnames(faustDF)
cellPops <- setdiff(cellPops,"0_0_0_0_0")
allFaustPops <- length(cellPops)
faustDF$parentCount <- apply(faustDF,1,sum)
faustDF$sampleName <- rownames(faustDF)
faustDF$subjectID <- rownames(faustDF)

#
#score FAUST and both flowSOM outputs as clusterings using F-measure (as well as other measures)
#
trueLabels <- Reduce(append,simmedExp[["labelsList"]])
trueNumericLabels <- as.numeric(as.factor(trueLabels))
#warnings suppressed because normalized mutual information and other information measures 
#do not work with bigInt.
flowsomScoresAll <- suppressWarnings(myExternalValidation(trueNumericLabels,
                                                          clusterLabels,
                                                          methodName=paste0("flowSOM_All_Iter",simIter)))
overpartitionedScoresAll <- suppressWarnings(myExternalValidation(trueNumericLabels,
                                                          opClusterLabels,
                                                          methodName=paste0("overpartitioned_All_Iter",simIter)))

faustSamplePaths <- list.files(paste0(jobPath,"/faustData/sampleData"),full.names=TRUE)
firstSample <- TRUE
for (samplePath in faustSamplePaths) {
    fa <- read.table(paste0(samplePath,"/faustAnnotation.csv"),
                     header=F,
                     sep="`")[,1]
    if (firstSample) {
        faustClustering <- as.character(fa)
        firstSample <- FALSE
    }
    else {
        faustClustering <- append(faustClustering,as.character(fa))
    }
}
faustScoresAll <- suppressWarnings(myExternalValidation(trueNumericLabels,
                                                        as.numeric(as.factor(faustClustering)),
                                                        methodName=paste0("FAUST_All_Iter_",simIter)))
#
#look up the subset of FAUST with phenotypes
#
faustPhenoLookup <- which(faustClustering != "0_0_0_0_0")
faustScoresPheno <- suppressWarnings(myExternalValidation(trueNumericLabels[faustPhenoLookup],
                                                          as.numeric(as.factor(faustClustering[faustPhenoLookup])),
                                                          methodName=paste0("FAUST_Pheno_Subset_Iter_",simIter)))

flowsomScoresPheno<- suppressWarnings(myExternalValidation(trueNumericLabels[faustPhenoLookup],
                                                           clusterLabels[faustPhenoLookup],
                                                           methodName=paste0("flowSOM_Pheno_Subset_Iter_",simIter)))

overpartitionedScoresPheno<- suppressWarnings(myExternalValidation(trueNumericLabels[faustPhenoLookup],
                                                                   opClusterLabels[faustPhenoLookup],
                                                           methodName=paste0("overpartitioned_Pheno_Subset_Iter_",simIter)))

clusterScores <- cbind(flowsomScoresAll,overpartitionedScoresAll,faustScoresAll,
                       flowsomScoresPheno,overpartitionedScoresPheno,faustScoresPheno)

#
#Derive a non-deterministic response based on the "spiked" population.
#spikedInPopStatusDF[,"spikedInPop"]==0 -> the population is not enriched.
#
spikedInPopStatusDF <- data.frame(spikedInPop=as.numeric(unlist(simmedExp[["spikedInPopList"]])),
                             sampleName=names(simmedExp[["spikedInPopList"]]),
                             stringsAsFactors=FALSE)
nonresp <- which(allSubjects %in% spikedInPopStatusDF[which(spikedInPopStatusDF[,"spikedInPop"]==0),"sampleName"])
#get medians for experiment to use for flowsom annotations
experimentMedians <- apply(exprsMat,2,median)
allFaustResults <- allFlowsomResults <- allOverpartitionedResults <- list()
for (responseRate in seq(0.5,0.8,by=0.05)) {
    print(paste0("Current response rate: ",responseRate))
    #
    #modify NITERS to change number of times response status is sampled.
    #
    modelIters <- NITERS <- 50
    flowsomMaxBestFDR <- flowsomMaxPropFDR <- flowsomMaxNameFDR <- numFlowsomPops <- bestFlowsomFDR <- c()
    overpartMaxBestFDR <- overpartMaxPropFDR <- overpartMaxNameFDR <-numOverpartitionedPops <- bestOverpartitionedFDR <- c()
    faustPvalList <- flowsomPvalList <- flowsomAnnotations <- list()
    overpartitionedPvalList <- overpartitionedAnnotations <- list()
    spikedFAUSTFDR <- numFAUSTPops <- bestFAUSTFDR <- c()
    faustAnnotations <- list()
    oracleBest <- oracleCount <- oracleProp<- list()
    opBest <- opCount <- opProp<- list()
    while (NITERS) {
        responseProb <- rep(responseRate,length(allSubjects))
        responseProb[nonresp] <- (1-responseRate)
        names(responseProb) <- allSubjects
        responderStatus <- sapply(responseProb,function(x){rbinom(1,1,x)})
        resDF <- data.frame(modelRT=as.factor(as.numeric(responderStatus)),
                            subjectID=names(responderStatus),
                            stringsAsFactors=FALSE) 
        #
        #join outcome to flowSOM and model counts
        #
        flowsomCountDF <- as.data.frame(flowsomCountMatrix)
        flowsomCellPops <- colnames(flowsomCountDF)
        flowsomCountDF$parentCount <- apply(flowsomCountDF,1,sum)
        flowsomCountDF$sampleName <- rownames(flowsomCountDF)
        flowsomCountDF$subjectID <- rownames(flowsomCountDF)
        flowsomModelDF <- left_join(flowsomCountDF,resDF,by="subjectID")
        flowsomPvalVec <- c()
        for (cellPop in flowsomCellPops) {
            mdf <- data.frame(resType=as.factor(flowsomModelDF[,"modelRT"]),
                              subjectID=as.factor(flowsomModelDF[,"subjectID"]),
                              parentCount=flowsomModelDF[,"parentCount"],
                              childCount=flowsomModelDF[,cellPop])
            suppressMessages(pv <- safeMod(mdf))
            if (!is.na(pv)) {
                flowsomPvalVec <- append(flowsomPvalVec,pv[1])
                names(flowsomPvalVec)[length(flowsomPvalVec)] <- paste0(cellPop)
            }
        }
        flowsomPvalList <- append(flowsomPvalList,list(flowsomPvalVec))
        flowsomFDR <- p.adjust(flowsomPvalVec,"BH")
        #
        #record flowsom cluster with containing absolute most of spiked in
        #
        if (any(modMaxCountOracleNames %in% names(flowsomFDR))) {
            flowsomMaxNameFDR <- append(flowsomMaxNameFDR,
                                        as.numeric(flowsomFDR[intersect(modMaxCountOracleNames,names(flowsomFDR))]))
            oracleCount <- append(oracleCount,1)
        }
        else {
            oracleCount <- append(oracleCount,0)
        }
        #
        #record flowsom cluster with containing proportionally most of spiked in
        #
        if (any(modMaxPropOracleNames %in% names(flowsomFDR))) {
            flowsomMaxPropFDR <- append(flowsomMaxPropFDR,
                                        as.numeric(flowsomFDR[intersect(modMaxPropOracleNames,names(flowsomFDR))]))
            oracleProp <- append(oracleProp,1)
        }
        else {
            oracleProp <- append(oracleProp,0)
        }
        #
        #record best of two
        #
        if (any(modMaxCountOracleNames %in% names(flowsomFDR)) || any(modMaxPropOracleNames %in% names(flowsomFDR))) {
            if (any(modMaxCountOracleNames %in% names(flowsomFDR)) && any(modMaxPropOracleNames %in% names(flowsomFDR))) {
                flowsomMaxBestFDR <- append(flowsomMaxBestFDR,
                                            min(as.numeric(flowsomFDR[intersect(modMaxPropOracleNames,names(flowsomFDR))]),
                                                as.numeric(flowsomFDR[intersect(modMaxCountOracleNames,names(flowsomFDR))])))
            }
            else if (any(modMaxPropOracleNames %in% names(flowsomFDR))) {
                flowsomMaxBestFDR <- append(flowsomMaxBestFDR,
                                            as.numeric(flowsomFDR[intersect(modMaxPropOracleNames,names(flowsomFDR))]))
            }
            else {
                flowsomMaxBestFDR <- append(flowsomMaxBestFDR,
                                            as.numeric(flowsomFDR[intersect(modMaxCountOracleNames,names(flowsomFDR))]))
            }
            oracleBest <- append(oracleBest,1)
        }
        else {
            oracleBest <- append(oracleBest,0)
        }
        bestFlowsomFDR <- append(bestFlowsomFDR,min(flowsomFDR))
        bestFlowsomPops <- flowsomFDR[which(flowsomFDR==min(flowsomFDR))]
        numFlowsomPops <- append(numFlowsomPops,length(bestFlowsomPops))
        minName <- names(bestFlowsomPops)[1]
        minLookup <- which(clusterLabels==minName)
        flowsomMedians <- apply(exprsMat[minLookup,],2,median)
        bestLabelExample <- paste0(paste0(names(which(flowsomMedians >= experimentMedians)),"+"),collapse="")
        flowsomAnnotations <- append(flowsomAnnotations,bestLabelExample)
        #
        #join outcome to overpartitioned flowSOM and model counts
        #
        overpartitionedCountDF <- as.data.frame(overpartitionedCountMatrix)
        overpartitionedCellPops <- colnames(overpartitionedCountDF)
        overpartitionedCountDF$parentCount <- apply(overpartitionedCountDF,1,sum)
        overpartitionedCountDF$sampleName <- rownames(overpartitionedCountDF)
        overpartitionedCountDF$subjectID <- rownames(overpartitionedCountDF)
        overpartitionedModelDF <- left_join(overpartitionedCountDF,resDF,by="subjectID")
        overpartitionedPvalVec <- c()
        for (cellPop in overpartitionedCellPops) {
            mdf <- data.frame(resType=as.factor(overpartitionedModelDF[,"modelRT"]),
                              subjectID=as.factor(overpartitionedModelDF[,"subjectID"]),
                              parentCount=overpartitionedModelDF[,"parentCount"],
                              childCount=overpartitionedModelDF[,cellPop])
            suppressMessages(pv <- safeMod(mdf))
            if (!is.na(pv)) {
                overpartitionedPvalVec <- append(overpartitionedPvalVec,pv[1])
                names(overpartitionedPvalVec)[length(overpartitionedPvalVec)] <- paste0(cellPop)
            }
        }
        overpartitionedPvalList <- append(overpartitionedPvalList,list(overpartitionedPvalVec))
        overpartitionedFDR <- p.adjust(overpartitionedPvalVec,"BH")
        #
        #record overpart cluster with containing absolute most of spiked in
        #
        if (any(modMaxCountOverpartNames %in% names(overpartitionedFDR))) {
            overpartMaxNameFDR <- append(overpartMaxNameFDR,
                                         as.numeric(overpartitionedFDR[intersect(modMaxCountOverpartNames,
                                                                                 names(overpartitionedFDR))]))
            opCount <- append(opCount,1)
        }
        else {
            opCount <- append(opCount,0)
        }
        #
        #record overpart cluster with containing proportionally most of spiked in
        #
        if (any(modMaxPropOverpartNames %in% names(overpartitionedFDR))) {
            overpartMaxPropFDR <- append(overpartMaxPropFDR,
                                         as.numeric(overpartitionedFDR[intersect(modMaxPropOverpartNames,
                                                                                 names(overpartitionedFDR))]))
            opProp <- append(opProp,1)
        }
        else {
            opProp <- append(opProp,0)
        }
        #
        #record best of two
        #
        if (any(modMaxCountOverpartNames %in% names(overpartitionedFDR)) ||
            any(modMaxPropOverpartNames %in% names(overpartitionedFDR))) {
            if (any(modMaxCountOverpartNames %in% names(overpartitionedFDR)) &&
                any(modMaxPropOverpartNames %in% names(overpartitionedFDR))) {
                overpartMaxBestFDR <- append(overpartMaxBestFDR,
                                             min(as.numeric(overpartitionedFDR[intersect(modMaxPropOverpartNames,
                                                                                         names(overpartitionedFDR))]),
                                                as.numeric(overpartitionedFDR[intersect(modMaxCountOverpartNames,
                                                                                        names(overpartitionedFDR))])))
            }
            else if (any(modMaxPropOverpartNames %in% names(overpartitionedFDR))) {
                overpartMaxBestFDR <- append(overpartMaxBestFDR,
                                             as.numeric(overpartitionedFDR[intersect(modMaxPropOverpartNames,
                                                                                     names(overpartitionedFDR))]))
            }
            else {
                overpartMaxBestFDR <- append(overpartMaxBestFDR,
                                             as.numeric(overpartitionedFDR[intersect(modMaxCountOverpartNames,
                                                                                     names(overpartitionedFDR))]))
            }
            opBest <- append(opBest,1)
        }
        else {
            opBest <- append(opBest,0)
        }
 


        bestOverpartitionedFDR <- append(bestOverpartitionedFDR,min(overpartitionedFDR))
        bestOverpartitionedPops <- overpartitionedFDR[which(overpartitionedFDR==min(overpartitionedFDR))]
        numOverpartitionedPops <- append(numOverpartitionedPops,length(bestOverpartitionedPops))
        opMinName <- names(bestOverpartitionedPops)[1]
        opMinLookup <- which(opClusterLabels==opMinName)
        overpartitionedMedians <- apply(exprsMat[opMinLookup,],2,median)
        bestOpLabelExample <- paste0(paste0(names(which(overpartitionedMedians >= experimentMedians)),"+"),collapse="")
        overpartitionedAnnotations <- append(overpartitionedAnnotations,bestOpLabelExample)

        
        #
        #join responder status onto faust counts, fit GLMMs, and record statistics
        #
        modelDF <- left_join(faustDF,resDF,by="subjectID")
        mePvalVec <- c()
        for (cellPop in cellPops) {
            mdf <- data.frame(resType=as.factor(modelDF[,"modelRT"]),
                              subjectID=as.factor(modelDF[,"subjectID"]),
                              parentCount=modelDF[,"parentCount"],
                              childCount=modelDF[,cellPop])
            suppressMessages(pv <- safeMod(mdf))
            if (!is.na(pv)) {
                mePvalVec <- append(mePvalVec,pv[1])
                names(mePvalVec)[length(mePvalVec)] <- paste0(cellPop)
            }
        }
        faustPvalList <- append(faustPvalList,list(mePvalVec))
        qp <- p.adjust(mePvalVec,"BH")
        bestFAUSTFDR <- append(bestFAUSTFDR,min(qp))
        bestCorrelate <- names(which(qp==min(qp)))
        numFAUSTPops <- append(numFAUSTPops,length(bestCorrelate))
        if (spikedInFound) {
            spikedFAUSTFDR<- append(spikedFAUSTFDR,qp[spikedInPop])
        }
        else {
            spikedFAUSTFDR <- append(spikedFAUSTFDR,1)
        }
        faustAnnotations <- append(faustAnnotations,
                                   length(intersect(spikedInPop,bestCorrelate)))
        NITERS <- NITERS - 1
        print(paste0("Regression iterations remaining: ",NITERS))
    }
    #
    #record flowsom results
    #
    flAnnSum <- table(unlist(flowsomAnnotations))
    flMaxName <- max(flAnnSum)
    flMaxNameEx <- names(flAnnSum[which(flAnnSum==flMaxName)[1]])
    flMaxNameEx

    flowsomResults <- list(
        responseRate=responseRate,
        medPopNumflowSOM=median(numFlowsomPops),
        medFDRflowSOM=median(bestFlowsomFDR),
        medOracleMaxNameFDR=median(flowsomMaxNameFDR),
        medOracleMaxPropFDR=median(flowsomMaxPropFDR),
        medOracleMaxBestFDR=median(flowsomMaxBestFDR),
        oracleCountFreq=sum(unlist(oracleCount)),
        oraclePropFreq=sum(unlist(oracleProp)),
        oracleBestFreq=sum(unlist(oracleBest)),
        flMaxName=flMaxName,
        flMaxNameEx=flMaxNameEx,
        modelIters=modelIters,
        flowsomPvalList=flowsomPvalList
    )
    allFlowsomResults <- append(allFlowsomResults,list(flowsomResults))
    #
    #record overpartitioned results
    #
    opAnnSum <- table(unlist(overpartitionedAnnotations))
    opMaxName <- max(opAnnSum)
    opMaxNameEx <- names(opAnnSum[which(opAnnSum==opMaxName)[1]])
    opMaxNameEx
    overpartitionedResults <- list(
        responseRate=responseRate,
        medPopNumoverpartitioned=median(numOverpartitionedPops),
        medFDRoverpartitioned=median(bestOverpartitionedFDR),
        medOPMaxNameFDR=median(overpartMaxNameFDR),
        medOpMaxPropFDR=median(overpartMaxPropFDR),
        medOpMaxBestFDR=median(overpartMaxBestFDR),
        opCountFreq=sum(unlist(opCount)),
        opPropFreq=sum(unlist(opProp)),
        opBestFreq=sum(unlist(opBest)),
        opMaxName=opMaxName,
        opMaxNameEx=opMaxNameEx,
        modelIters=modelIters,
        overpartitionedPvalList=overpartitionedPvalList
    )
    allOverpartitionedResults <- append(allOverpartitionedResults,list(overpartitionedResults))
    #
    #record faust results
    #
    faustAnnSum <- sum(unlist(faustAnnotations))
    faustResults <- list(
        responseRate=responseRate,
        medianSpikedFDR=median(spikedFAUSTFDR),
        medPopNumFAUST=median(numFAUSTPops),
        medFDRFAUST=median(bestFAUSTFDR),
        faustMaxName=faustAnnSum,
        modelIters=modelIters,
        faustPvalList=faustPvalList
    )
    allFaustResults <- append(allFaustResults,list(faustResults))
}
#
#record all iteration results
#
simulationResults <- list(
    spikedInPop=spikedInPop,
    medianFaustPctUnclassified=medianFaustPctUnclassified,
    clusterScores=clusterScores,
    faustVarSelectionCount=faustVarSelectionCount,
    allFaustPops=allFaustPops,
    spikedCorrelation=spikedCorrelation,
    medFaustCorrelation=medFaustCorrelation,
    spikedInFound=spikedInFound,
    truePopsFound=truePopsFound,
    faustModeling=allFaustResults,
    flowsomModeling=allFlowsomResults,
    overpartitionedModeling=allOverpartitionedResults    
)
allSimulationResults <- append(allSimulationResults,list(simulationResults))
#
#decrement loop, clean up, and simulate another experiment unless termination signal found.
#
simIter <- simIter-1
unlink(faustDataPath,recursive=TRUE,force=FALSE)
#save intermediate data to monitor job progress.
saveRDS(allSimulationResults,tmpResultsPath)
if (file.exists(fileAbortPath)) {
    print("Job abort signal detected.")
    print(simIter)
    simIter <- 0
}
}
#
#simulation is complete. save final results.
#
saveRDS(allSimulationResults,resultsPath)
                                 
