#' Full Annotation Using Shape-constrained Trees
#'
#' makeAnnotationEmbedding produces an annotation embedding for samples analyzed
#' by FAUST.
#'
#' The annotation embedding is constructed by first windsorizing each sample
#' to the 1st and 99th percentile. Then, the observed measurements for each FAUST annotation
#' are normalized to mean 0 unit SD (if >1). Finally the normalized annotation measurements
#' are translated by coordinate onto annotation landmarks. 
#'
#' After transforming each sample in this way, all listed samples are concatenated.
#' The concatendated array is embedded using the UMAP algorithm with the qualitative
#' parameter settings suggested in doi:10.1038/nbt.4314.
#' 
#' @param projectPath A path to a directory containing a completed FAUST analysis. 
#' The projectPath directory is expected to have a faustData sub-directory, which
#' contains a completed faust analysis. 
#'
#' @param sampleNameVec A character vector listing samples to embed. Each entry
#' in the vector must match a directory name in projectPath/faustData/sampleData.
#'
#' @param vizType A string listing the type of annotation embedding to perform.
#' The default value is "uniform", which causes the distance between any two
#' annotations is uniform across coordinates. The other supported value is
#' "depthscore", which scales annotation coordinates by their normalized depth score.
#' 
#' @param embeddingSelectionQuantile A numeric value between 0 and 1. Indicates which
#' the quantile of the depth-score distribution to use for marker inclusion in the 
#' annoation embedding.
#'
#' @param embeddingDepthScoreThreshold A numeric value between 0 and 1. Indicates which
#' the depth-score threshold value to use for marker inclusion in the 
#' annoation embedding.
#' 
#' @return A data frame with the following columns:
#'
#' Column `umapX` and `umapY`: the x-y coordinates of the embedding.
#'
#' Column `sampleOfOrigin`: character vector listing the sample of origin.
#'
#' Column `faustLabels`: the faust annotation of the observation.
#'
#' Column `marker name`: the observed expression data, windsorized at the
#' 1st and 99th percentile, and scaled to the unit interval. 
#'
#' Column `marker name_annot`: the FAUST annotation for the marker.
#' 
#' @export
#' @md
#' @importFrom uwot umap
makeAnnotationEmbedding <- function(
                                    projectPath,
                                    sampleNameVec,
                                    vizType="uniform",
                                    embeddingSelectionQuantile=0.75,
                                    embeddingDepthScoreThreshold=0.01
                                    )
{
    require(uwot)
    if (!dir.exists(file.path(projectPath,"faustData"))) {
        print(paste0("No faustData detected at projectPath location: ",projectPath))
        stop("Update path or run FAUST before calling this function.")
    }
    if (!dir.exists(file.path(projectPath,"faustData","gateData"))) {
        print(paste0("No faustData/gateData detected at projectPath location: ",projectPath))
        stop("Update path or run FAUST before calling this function.")
    }
    if (!dir.exists(file.path(projectPath,"faustData","metaData"))) {
        print(paste0("No faustData/metaData detected at projectPath location: ",projectPath))
        stop("Update path or run FAUST before calling this function.")
    }
    if (!dir.exists(file.path(projectPath,"faustData","sampleData"))) {
        print(paste0("No faustData/sampleData detected at projectPath location: ",projectPath))
        stop("Update path or run FAUST before calling this function.")
    }
    allSampleFiles <- list.files(file.path(projectPath,"faustData","sampleData"))
    if (length(setdiff(sampleNameVec,allSampleFiles))) {
        print("The following samples are not detected in projectPath/faustData/sampleData")
        print(setdiff(sampleNameVec,allSampleFiles))
        stop("Update path before calling this function.")
    }
    if (length(intersect(sampleNameVec,allSampleFiles)) != length(sampleNameVec)) {
        print("Attempting to include a sample in the embedding that's not in projectPath/faustData/sampleData.")
        stop("Update path before calling this function.")
    }
    if (!(vizType %in% (c("uniform","depthscore")))) {
        stop("Unsupported vizType detected. Set vizType parameter to uniform or depthscore")
    }
    preppedSampleData <- .prepareSamplesForEmbed(
        sampleNameVec=sampleNameVec,
        projectPath=projectPath,
        visualizationType=vizType
    )
    preppedDataForEmbedding <- preppedSampleData$preparedDataset
    #gatingFiles <- list.files(file.path(projectPath,"faustData","gateData"))
    #selectedMarkersFN <- gatingFiles[grepl("selectedChannels",gatingFiles)]
    #selectedMarkers <- readRDS(file.path(projectPath,"faustData","gateData",selectedMarkersFN))
    allSelectedMarkers <- colnames(preppedDataForEmbedding)
    faustDepthScoreMat <- readRDS(file.path(projectPath,"faustData","metaData","scoreMat.rds"))
    embeddingScores <- apply(faustDepthScoreMat,2,function(x){quantile(x,probs=embeddingSelectionQuantile)})
    embeddingSelMarkers <- names(which(embeddingScores >= embeddingDepthScoreThreshold))
    embeddingMarkers <- intersect(embeddingSelMarkers,allSelectedMarkers)
    finalMatrixForEmbedding <- preppedDataForEmbedding[,embeddingMarkers,drop=FALSE]
    #
    #embed the prepared dataset.
    #
    #annoEmbed <- uwot::umap(preppedDataForEmbedding,n_neighbors=15,metric="euclidean",min_dist=0.2,verbose=FALSE)
    annoEmbed <- uwot::umap(finalMatrixForEmbedding,n_neighbors=15,metric="euclidean",min_dist=0.2,verbose=FALSE)
    umapData <- data.frame(
        umapX=annoEmbed[,1],
        umapY=annoEmbed[,2],
        stringsAsFactors=FALSE
    )
    #
    #attach sample lookups and selected faust clusters
    #
    umapData$sampleOfOrigin <- as.character(preppedSampleData$sampleLookups)
    umapData$faustLabels <- preppedSampleData$clusterLabels
    #
    #attach the observed expression
    #
    rawExprsMat <- preppedSampleData$rawExpressionMatrix
    umapDataOutPrep <- cbind(umapData,rawExprsMat)
    #
    #attach the windsorized expression
    #
    windExpressionData <- preppedSampleData$windsorizedExpressionMatrix
    umapDataOutPrep2 <- cbind(umapDataOutPrep,windExpressionData)
    #
    #attach the exhaustive faust annotations
    #
    exAnnot <- preppedSampleData$exhaustiveAnnotations
    umapDataOut <- cbind(umapDataOutPrep2,exAnnot)
    return(umapDataOut)
}


.windsorizeAndScaleExpression <- function(eVec) {
    windsorLow <- as.numeric(quantile(eVec,probs=c(0.01)))
    windsorHigh <- as.numeric(quantile(eVec,probs=c(0.99)))
    windVec <- eVec
    windVec[which(windVec <= windsorLow)] <- windsorLow
    windVec[which(windVec >= windsorHigh)] <- windsorHigh
    windVecOut <- ((windVec-min(windVec))/(max(windVec)-min(windVec)))
    return(windVecOut)
}

.windsorizeExpressionForEmbed <- function(eVec) {
    windsorLow <- as.numeric(quantile(eVec,probs=c(0.01)))
    windsorHigh <- as.numeric(quantile(eVec,probs=c(0.99)))
    windVec <- eVec

    lowLookup <- which(windVec <= windsorLow)
    lowRep <- rep(windsorLow,length(lowLookup))
    #
    #add perturbations to avoid exact ties across very rare
    #annotation groups
    #
    lowRep <- lowRep + rnorm(n=length(lowLookup),mean=0,sd=0.01)
    windVec[lowLookup] <- lowRep 

    highLookup <- which(windVec >= windsorHigh)
    highRep <- rep(windsorHigh,length(highLookup))
    highRep <- highRep + rnorm(n=length(highLookup),mean=0,sd=0.01)
    windVec[highLookup] <- highRep

    return(windVec)
}

.standardizeVectorForAnnotScale <- function(vec) {
    vecMean <- mean(vec)
    vecSD <- sd(vec)
    if (vecSD == 0) {
        vecOut <- (vec-vecMean)
    }
    else {
        vecOut <- ((vec-vecMean)/vecSD)
    }
        
    return(vecOut)
}

.createListForAnnotRecords <- function(alList) 
{
    annotRecordDim <- 2*sum(unlist(lapply(alList,length)))
    annotRecords <- vector(mode="list",annotRecordDim)
    recordNameTemplate <- unlist(lapply(names(alList),function(x){paste0(x,"%",alList[[x]])}))
    names(annotRecords) <- c(paste0(recordNameTemplate,"%Low"),paste0(recordNameTemplate,"%High"))
    for (arName in names(annotRecords)) {
        annotRecords[[arName]] <- c(NA)
    }
    return(annotRecords)
}

.getTargetThresholds <- function(selectedChannels,alList)
{
    targetThresholdList <- list()
    for (selMarker in selectedChannels) {
        scaledLevels <- sort(alList[[selMarker]])
        markerThresholdVec <- c()
        for (tNum in seq(1,(length(scaledLevels)-1))) {
            markerThresholdVec <- append(markerThresholdVec,
                                         mean(scaledLevels[c(tNum,(tNum+1))]))
        }
        targetThresholdList <- append(targetThresholdList,list(markerThresholdVec))
        names(targetThresholdList)[length(targetThresholdList)] <- selMarker
    }
    return(targetThresholdList)
}
        
.getQVforAnnotEmbed <- function(nameSource,sampleAnnotRecords,markerQuantilesSource) {
    markerName <- strsplit(nameSource,"%")[[1]][1]
    markerQV <- markerQuantilesSource[markerName]
    return(as.numeric(quantile(sampleAnnotRecords[[nameSource]],probs=markerQV,na.rm=TRUE)))
}

.computeScalingFactors <- function(targetThreshList,alList,minScalingValues,maxScalingValues) {
    finalScalingFactors <- list()
    for (markerName in names(targetThreshList)) {
        currentThresholds <- targetThreshList[[markerName]]
        currentCenters <- alList[[markerName]]
        currentScenarios <- paste0(markerName,"%",seq(length(currentCenters)))
        for (scenarioNum in seq(length(currentScenarios))) {
            scenarioDesc <- currentScenarios[scenarioNum]
            if (scenarioNum==1) {
                #the lowest expression group
                threshHigh <- currentThresholds[scenarioNum]
                translate <- currentCenters[scenarioNum]
                scalingValHigh <- maxScalingValues[[scenarioDesc]]
                scenarioOutcome <- ((threshHigh-translate)/scalingValHigh)
            }
            else if (scenarioNum == length(currentScenarios)) {
                #the highest expression group
                threshLow <- currentThresholds[(scenarioNum-1)]
                translate <- currentCenters[scenarioNum]
                scalingValLow <- minScalingValues[[scenarioDesc]]
                scenarioOutcome <- ((threshLow-translate)/scalingValLow)
            }
            else {
                #an intermediate expression group. have to handle both scenarios.
                translate <- currentCenters[scenarioNum]

                threshHigh <- currentThresholds[scenarioNum]
                threshLow <- currentThresholds[(scenarioNum-1)]

                scalingValHigh <- maxScalingValues[[scenarioDesc]]
                scalingValLow <- minScalingValues[[scenarioDesc]]

                highOutcome <- ((threshHigh-translate)/scalingValHigh)
                lowOutcome <- ((threshLow-translate)/scalingValLow)
                scenarioOutcome <- min(highOutcome,lowOutcome)
            }
            finalScalingFactors <- append(finalScalingFactors,list(scenarioOutcome))
            names(finalScalingFactors)[length(finalScalingFactors)] <- scenarioDesc
        }
    }
    return(finalScalingFactors)
}

.prepareSampleForEmbed <- function(
                                  projectPath,
                                  sampleName,
                                  alList,
                                  localAnnotRecords,
                                  targetThresholdList,
                                  selectedChannels,
                                  markerQuantilesLow,
                                  markerQuantilesHigh
                                  )
{
    #
    #load the observed expression matrix
    #
    exprsMatIn <- readRDS(file.path(projectPath,"faustData","sampleData",sampleName,"exprsMat.rds"))
    exprsMat <- exprsMatIn[,selectedChannels,drop=FALSE]
    #
    #load the per-cell annotation matrix
    #
    amIn <- read.table(file.path(projectPath,"faustData","sampleData",sampleName,"annotationMatrix.csv"),
                       header=FALSE,
                       sep=",",
                       stringsAsFactors=FALSE)
    colnames(amIn) <- selectedChannels
    fullAnnVec <- apply(amIn,1,function(x){paste0(x,collapse="~")})
    uniqFullAnnVec <- unique(fullAnnVec)
    #
    #create a modified expression matrix scaled by faust cluster
    #
    clusterScaleExprsMat <- apply(exprsMat[,selectedChannels,drop=FALSE],2,.windsorizeExpressionForEmbed)
    for (ufa in uniqFullAnnVec) {
        ufaLookup <- which(fullAnnVec==ufa)
        clusterExpression <- clusterScaleExprsMat[ufaLookup,,drop=FALSE]
        if (nrow(clusterExpression) == 1) {
            scaledClusterExprs <- clusterExpression
            #set to a random deviate centered at 0. updated from setting to exact 0
            #to avoid scenario in rare samples where every cluster associated with 
            #an annotation grouping consists of a single observation
            scaledClusterExprs[1,] <- rnorm(ncol(clusterExpression),sd=1e-7)
        }
        else {
            scaledClusterExprs <- apply(clusterExpression,2,.standardizeVectorForAnnotScale)
        }
        clusterMins <- apply(scaledClusterExprs,2,min)
        clusterMaxs <- apply(scaledClusterExprs,2,max)
        annotType <- paste0(selectedChannels,"%",strsplit(ufa,"~")[[1]])
        #
        #update cluster statistics
        #
        for (atNum in seq(length(annotType))) {
            curAnnotType <- annotType[[atNum]]
            localAnnotRecords[[paste0(curAnnotType,"%Low")]] <- append(localAnnotRecords[[paste0(curAnnotType,"%Low")]],clusterMins[atNum])
            localAnnotRecords[[paste0(curAnnotType,"%High")]] <- append(localAnnotRecords[[paste0(curAnnotType,"%High")]],clusterMaxs[atNum])
        }
        #
        #update cell-level expression matrix with scaled values
        #
        for (obsNum in seq(length(ufaLookup))) {
            updateNum <- ufaLookup[obsNum]
            clusterScaleExprsMat[updateNum,] <- scaledClusterExprs[obsNum,]
        }
    }
    lowNames <- names(localAnnotRecords)[grepl("Low",names(localAnnotRecords))]
    highNames <- names(localAnnotRecords)[grepl("High",names(localAnnotRecords))]
    minScalingValues <- lapply(lowNames,function(x){.getQVforAnnotEmbed(x,localAnnotRecords,markerQuantilesLow)})
    names(minScalingValues) <- gsub("%Low","",lowNames)
    maxScalingValues <- lapply(highNames,function(x){.getQVforAnnotEmbed(x,localAnnotRecords,markerQuantilesHigh)})
    names(maxScalingValues) <- gsub("%High","",highNames)
    #
    #compute the scaling for each category
    #
    finalScalingFactors <- .computeScalingFactors(
        targetThreshList=targetThresholdList,
        alList=alList,
        minScalingValues=minScalingValues,
        maxScalingValues=maxScalingValues
    ) 
    #
    #finally update the cluster expression with the scalinings, and translate the clusters.
    #
    for (ufa in uniqFullAnnVec) {
        ufaLookup <- which(fullAnnVec==ufa)
        clusterExpression <- clusterScaleExprsMat[ufaLookup,,drop=FALSE]
        annotType <- paste0(selectedChannels,"%",strsplit(ufa,"~")[[1]])
        #
        #scale cluster by final value and translate it to depth score scale.
        #
        for (atNum in seq(length(annotType))) {
            curAnnotType <- annotType[[atNum]]
            #
            #this is a bug if the marker name uses the character '%'
            #
            scenarioData <-strsplit(curAnnotType,"%")[[1]]
            markerName <- scenarioData[[1]]
            currentCenters <- alList[[markerName]]
            translate <- currentCenters[as.numeric(scenarioData[[2]])]
            clusterScaling <- finalScalingFactors[[curAnnotType]]
            perturbations <- rnorm(nrow(clusterExpression),sd=1e-4)
            clusterExpression[,atNum] <- ((clusterScaling*clusterExpression[,atNum]) + perturbations + translate)
        }
        #
        #update cell-level expression matrix with scaled values
        #
        for (obsNum in seq(length(ufaLookup))) {
            updateNum <- ufaLookup[obsNum]
            clusterScaleExprsMat[updateNum,] <- clusterExpression[obsNum,]
        }
    }
    return(clusterScaleExprsMat)
}

.prepareSamplesForEmbed <- function(sampleNameVec,projectPath,visualizationType) {
    #
    #load the sample expression and annotations
    #
    gatingFiles <- list.files(file.path(projectPath,"faustData","gateData"))
    selCFN <- gatingFiles[grepl("_selectedChannels.rds",gatingFiles)]
    selC <- readRDS(file.path(projectPath,"faustData","gateData",selCFN))
    #
    #bring in the annotation map and create storage list
    #
    cnMap <- readRDS(file.path(projectPath,"faustData","metaData","colNameMap.rds"))
    cnTemplate <- cnMap[1,"faustColNames"]
    cnData <- strsplit(cnTemplate,"~")[[1]]
    possibleAnnotationLevels <- as.numeric(cnData[seq(from=3,to=length(cnData),by=3)])
    channelLabels <- cnData[seq(from=1,to=length(cnData),by=3)]
    if (!(TRUE==identical(channelLabels,selC))) {
        stop("Error. Generated an annotation label for an unselected marker")
    }
    annotLevelList <- lapply(possibleAnnotationLevels,function(x){seq(x)})
    names(annotLevelList) <- channelLabels
    annotRecords <- .createListForAnnotRecords(alList=annotLevelList)
    #
    #scale annotations relative to depth scores
    #
    scoreMat <- readRDS(file.path(projectPath,"faustData","metaData","scoreMat.rds"))
    rawScoreMultipliers <- apply(scoreMat[,selC,drop=FALSE],2,sum)
    scoreMultipliers <- rawScoreMultipliers/nrow(scoreMat)
    #normalize relative to min(AUC,0.1) and scale by golden ratio (1.618)
    scoreMultipliers <- (1.618*(scoreMultipliers/min(min(scoreMultipliers),0.1)))
    for (selMarker in names(scoreMultipliers)) {
        scaledLevels <- sort(annotLevelList[[selMarker]])
        if (visualizationType=="uniform") {
            annotLevelList[[selMarker]] <- (scaledLevels * 10)
        }
        else {
            #depth score scaling
            annotLevelList[[selMarker]] <- (scaledLevels * scoreMultipliers[[selMarker]])
        }
    }
    quantileScoreWeights <- rawScoreMultipliers/sum(rawScoreMultipliers)
    quantileScoreWeights <- c(0,quantileScoreWeights[-length(quantileScoreWeights)])
    names(quantileScoreWeights) <- names(rawScoreMultipliers)
    markerQuantilesLow <- cumsum(quantileScoreWeights)*0.025
    markerQuantilesHigh <- 1-markerQuantilesLow
    #
    #get the target threshold(s) for each marker
    #
    targetThresholdList <- .getTargetThresholds(selectedChannels=selC,alList=annotLevelList)
    firstSample <- TRUE
    for(sName in sampleNameVec) {
        preparedSample <- .prepareSampleForEmbed(
            projectPath=projectPath,
            sampleName=sName,
            alList=annotLevelList,
            localAnnotRecords=annotRecords,
            targetThresholdList=targetThresholdList,
            selectedChannels=selC,
            markerQuantilesLow=markerQuantilesLow,
            markerQuantilesHigh=markerQuantilesHigh
        )
        #
        #again load the expression data for reporting.
        #
        exprsMatIn <- readRDS(file.path(projectPath,"faustData","sampleData",sName,"exprsMat.rds"))
        #
        #windsorize and scale the selected markers
        #
        exprsMat <- exprsMatIn[,selC,drop=FALSE]
        wExprsMat <- apply(exprsMat,2,.windsorizeAndScaleExpression)
        windColNames <- paste0(colnames(wExprsMat),"_Windsorized")
        colnames(wExprsMat) <- windColNames
        #
        #similarly, again load the per-cell annotation matrix
        #
        amIn <- read.table(file.path(projectPath,"faustData","sampleData",sName,"annotationMatrix.csv"),
                           header=FALSE,
                           sep=",",
                           stringsAsFactors=FALSE)
        colnames(amIn) <- selC
        #apply interpretable labels to encodings
        amUpdate <- apply(amIn,2,as.character)
        for (cmn in selC) {
            cNumLevels <- length(annotLevelList[[cmn]])
            cUpdateVec <- amUpdate[,cmn,drop=TRUE]
            if (cNumLevels == 2) {
                cUpdateVec <- gsub("1","-",cUpdateVec)
                cUpdateVec <- gsub("2","+",cUpdateVec)
            }
            else if (cNumLevels == 3) {
                cUpdateVec <- gsub("1","-",cUpdateVec)
                cUpdateVec <- gsub("2","Dim",cUpdateVec)
                cUpdateVec <- gsub("3","Bright",cUpdateVec)
            }
            else if (cNumLevels == 4) {
                cUpdateVec <- gsub("1","-",cUpdateVec)
                cUpdateVec <- gsub("2","MedLow",cUpdateVec)
                cUpdateVec <- gsub("3","MedHigh",cUpdateVec)
                cUpdateVec <- gsub("4","VeryBright",cUpdateVec)
            }
            else {
                print(paste0("Greater than 4 expression categories detected for ",cmn))
            }
            amUpdate[,cmn] <- cUpdateVec
        }
        namColNames <- paste0(colnames(amUpdate),"_faust_annotation")
        colnames(amUpdate) <- namColNames
        #
        #also load the faust cluster labels for reporting.
        #
        fa <- read.table(file.path(projectPath,"faustData","sampleData",sName,"faustAnnotation.csv"),
                         header=FALSE,
                         sep="`",
                         stringsAsFactors=FALSE)[,1]
        fa <- gsub("~1~2~","-",fa)
        fa <- gsub("~2~2~","+",fa)
        fa <- gsub("~1~3~","-",fa)
        fa <- gsub("~2~3~","Dim",fa)
        fa <- gsub("~3~3~","Bright",fa)
        fa <- gsub("~1~4~","-",fa)
        fa <- gsub("~2~4~","MedLow",fa)
        fa <- gsub("~3~4~","MedHigh",fa)
        fa <- gsub("~4~4~","VeryBright",fa)
        #make a record of sample provenance
        sLookupVec <- rep(sName,nrow(preparedSample))
        if (firstSample) {
            preparedData <- preparedSample
            allSampleLookups <- sLookupVec
            allExprsData <- exprsMatIn
            allWindData <- wExprsMat
            allFA <- fa
            allExhaustiveAnn <- amUpdate
            firstSample <- FALSE
        }
        else {
            preparedData <- rbind(preparedData,preparedSample)
            allSampleLookups <- append(allSampleLookups,sLookupVec)
            allExprsData <- rbind(allExprsData,exprsMatIn)
            allWindData <- rbind(allWindData,wExprsMat)
            allFA <- append(allFA,fa)
            allExhaustiveAnn <- rbind(allExhaustiveAnn,amUpdate)
        }
    }
    outList <- list(
        `preparedDataset`= preparedData,
        `sampleLookups`= allSampleLookups,
        `rawExpressionMatrix`=allExprsData,
        `windsorizedExpressionMatrix`=allWindData,
        `clusterLabels`=allFA,
        `exhaustiveAnnotations`=allExhaustiveAnn 
    )
    return(outList)
}
