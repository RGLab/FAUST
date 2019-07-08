.getEmpiricalDepthCounts <- function(maxDepth,annotationForest)
{
    #count the number of times a marker is cut at each possible depth
    depthCounts <- rep(0,maxDepth)
    names(depthCounts) <- seq(maxDepth)
    for (channelNum in seq(1,length(annotationForest),by=5)) {
        channelName <- names(annotationForest)[[channelNum]]
        #trim strings that can be appended to channel name in the c++.
        channelName <- gsub("_cutPoints","",channelName)
        channelName <- gsub("_numCuts","",channelName)
        channelName <- gsub("_cutDepth","",channelName)
        channelName <- gsub("_nodePathScore","",channelName)
        channelName <- gsub("_nodePopSize","",channelName)
        channelCuts <- annotationForest[[paste0(channelName,"_numCuts")]]
        uniqCuts <- sort(unique(channelCuts))
        if (length(uniqCuts)>0) {
            channelDepths <- annotationForest[[paste0(channelName,"_cutDepth")]]
            for (cutNum in uniqCuts) {
                #the c convention: depth number is repeated in the number of 
                #cutpoints times in the annotation forest.
                #need to lookup up the subCuts to prevent over counting
                cutLookups <- which(channelCuts == cutNum)
                subCuts <- cutLookups[seq(1,length(cutLookups),by=cutNum)]
                currentDepths <- channelDepths[subCuts]
                for (cDepthNum in seq(maxDepth)) {
                    depthCounts[[cDepthNum]] <-  (depthCounts[[cDepthNum]] + length(which(currentDepths == cDepthNum)))
                }
            }
        }
    }
    return(depthCounts)
}

.getScoreListForChannel <- function(channel,annotationForest,eChannelSize,depthNorm,uniqNumCuts)
{
    nodeScoreList <- list()
    for (cpNum in uniqNumCuts) {
        #in addition to the empirical depth penalty,                                                                                                                                                            
        #a node is penalized by the cumulative product of (1-dipTestPvalue)                                                                                                                                            
        #along the path that led there. it also is penalized by the proportion                                                                                                                                         
        #of the nodes parent population to the root population in the forest.                                                                                                                                          
        allIndex <- which(annotationForest[[paste0(channel,"_numCuts")]]==cpNum)
        weightLookups <- allIndex[seq(1,length(allIndex),by=cpNum)]
        nodeSizes <- annotationForest[[paste0(channel,"_nodePopSize")]][weightLookups]
        popWeights <- nodeSizes/eChannelSize
        nodeDepths <- annotationForest[[paste0(channel,"_cutDepth")]][weightLookups]
        depthWeights <- depthNorm[nodeDepths]
        nodePathScores <- annotationForest[[paste0(channel,"_nodePathScore")]][weightLookups]
        nodeScores <- nodePathScores*popWeights*depthWeights
        nodeScoreList <- append(nodeScoreList,list(nodeScores))
        names(nodeScoreList)[length(nodeScoreList)] <- cpNum
    }
    return(nodeScoreList)
}


.parseAnnotationForest <- function(annotationForest,effectiveChannelSizes) {
    channelStrPrep <- names(annotationForest)[seq(from=1,to=length(annotationForest),by=5)]
    allChannels <- as.character(sapply(channelStrPrep,function(x){gsub("_cutPoints","",x)}))
    channelDataList <- list()
    #if a channel is never cut, slot assigned -Inf                                                                                                                                                                             
    channelDepths <- lapply(allChannels,function(x){sort(unique(annotationForest[[paste0(x,"_cutDepth")]]))})
    maxDepth <- suppressWarnings(max(unlist(lapply(channelDepths,function(x){max(x)}))))
    if (maxDepth < 0) {
        #the annotationForest is empty. return a null parse.
        for (channel in allChannels) {
            channelData <- list(NA,NA)
            names(channelData) <- c("channelScore","gates")
            channelDataList <- append(channelDataList,list(channelData))
            names(channelDataList)[length(channelDataList)] <- channel
        }
        return(channelDataList) 
    }
    depthCounts <- .getEmpiricalDepthCounts(
        maxDepth=maxDepth,
        annotationForest=annotationForest
    )
    #each depth of the tree is normalized to sum to one with the empical counts                                                     
    depthNorm <- 1/depthCounts
    for (channel in allChannels) {
        #eChannelSize is the number of active cells in a channel that were used to estimate the density.
        eChannelSize <- effectiveChannelSizes[[channel]]
        uniqNumCuts <- unique(annotationForest[[paste0(channel,"_numCuts")]])
        if (length(uniqNumCuts) > 0) {
            nodeScoreList <- .getScoreListForChannel(
                channel=channel,
                annotationForest=annotationForest,
                eChannelSize=eChannelSize,
                depthNorm=depthNorm,
                uniqNumCuts=uniqNumCuts
            )
            cutScores <- unlist(lapply(nodeScoreList,sum))
            channelScore <- sum(cutScores)
            if (channelScore > 0) {
                bestCutPointScore <- max(cutScores)
                cutPointSelection <- names(which(cutScores==bestCutPointScore))[1] #if tied channels, pick min
                nodeWeights <- nodeScoreList[[cutPointSelection]]/sum(nodeScoreList[[cutPointSelection]])
                selNumCuts <- as.numeric(cutPointSelection)
                selIndex <- which(annotationForest[[paste0(channel,"_numCuts")]]==selNumCuts)
                gates <- c()
                for (cutNum in seq(as.numeric(cutPointSelection))) {
                    gateLookups <- selIndex[seq(cutNum,length(selIndex),by=selNumCuts)]
                    gateLocations <- annotationForest[[paste0(channel,"_cutPoints")]][gateLookups]
                    gateLocation <- stats::weighted.mean(gateLocations,nodeWeights)
                    gates <- append(gates,gateLocation)
                }
                channelData <- list(channelScore,gates)
            }
            else {
                channelData <- list(0,NA)
            }
        }
        else {
            channelData <- list(0,NA)
        }
        names(channelData) <- c("channelScore","gates")
        channelDataList <- append(channelDataList,list(channelData))
        names(channelDataList)[length(channelDataList)] <- channel
    }
    #normalize scores to one
    maxChannelScore <- max(unlist(lapply(channelDataList,function(x){x$channelScore})))
    for (i in seq(length(channelDataList))) {
        channelDataList[[i]]$channelScore <- ((channelDataList[[i]]$channelScore)/maxChannelScore)
    }
    return(channelDataList)
}

