.selectChannels <- function(parentNode,
                            analysisMap,
                            depthScoreThreshold,
                            selectionQuantile,
                            projectPath) {
    uniqueLevels <- unique(analysisMap[,"analysisLevel"])
    firstForest <- TRUE
    for (aLevel in uniqueLevels) {
        if (file.exists(file.path(normalizePath(projectPath),
                                  "faustData",
                                  "levelData",
                                  aLevel,
                                  paste0(parentNode,"_pAnnF.rds"))))
        {
            afIn <- readRDS(file.path(normalizePath(projectPath),
                                      "faustData",
                                      "levelData",
                                      aLevel,
                                      paste0(parentNode,"_pAnnF.rds")))
            forestScores <- unlist(lapply(afIn,function(x){x$channelScore}))
            rafIn <- readRDS(file.path(normalizePath(projectPath),
                                       "faustData",
                                       "levelData",
                                       aLevel,
                                       paste0(parentNode,"_annF.rds")))
            gd <- rafIn[["gateData"]]
            firstDepth <- unlist(lapply(seq(3,length(gd),by=5),
                                        function(x){ifelse(length(gd[[x]])==0,Inf,min(gd[[x]]))}))
            names(firstDepth) <- names(forestScores)
            if (firstForest) {
                depthMat <- scoreMat <- matrix(nrow=0,ncol=length(forestScores))
                colnames(scoreMat) <- names(forestScores)
                colnames(depthMat) <- names(forestScores)
                firstForest <- FALSE
            }
            scoreMat <- rbind(scoreMat,forestScores)
            rownames(scoreMat)[length(rownames(scoreMat))] <- aLevel
            depthMat <- rbind(depthMat,firstDepth)
            rownames(depthMat)[length(rownames(depthMat))] <- aLevel
        }
        else {
            print(paste0("No parsed forest detected at aLevel: ",aLevel))
        }
    }
    if (firstForest) {
        return(c())
    }
    else {
        for (column in colnames(scoreMat)) {
            naLookup <- which(is.na(scoreMat[,column]))
            if (length(naLookup)) {
                scoreMat[naLookup,column] <- 0
            }
        }
        saveRDS(scoreMat,
                file.path(normalizePath(projectPath),
                          "faustData",
                          "metaData",
                          "scoreMat.rds"))
        saveRDS(depthMat,
                file.path(normalizePath(projectPath),
                          "faustData",
                          "metaData",
                          "depthMat.rds"))
        quantileScore <- apply(scoreMat,2,function(x){as.numeric(quantile(x,probs=c(selectionQuantile)))})
        scoreNames <- names(which(quantileScore >= depthScoreThreshold))
        quantileDepth <- apply(depthMat,2,function(x){as.numeric(quantile(x,probs=c((1-selectionQuantile))))})
        depthNames <- names(which(quantileDepth <= 2))
        selectedNames <- scoreNames 
        #outScore <- quantileScore[selectedNames]
        #
        #Order selected markers relative to their total depth score across the experiment.
        #
        outScore <- apply(scoreMat[,selectedNames,drop=FALSE],2,sum)
        outNames <- names(sort(outScore,decreasing=TRUE))
        return(outNames)
    }
}

