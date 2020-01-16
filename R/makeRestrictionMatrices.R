.makeRestrictionMatrices <- function(samplesInExp,
                                     analysisMap,
                                     channelBounds="",
                                     projectPath=".",
                                     debugFlag
                                     ) {
    if (!file.exists(file.path(normalizePath(projectPath),
                               "faustData",
                               "metaData",
                               "madeResMats.rds"))) {
        uniqueHierarchyNames <- names(table(analysisMap[,"impH",drop=TRUE]))
        uniqueHNum <- length(uniqueHierarchyNames)
        channelBoundsUsedByFAUST <- channelBounds
        if (channelBounds=="") {
            #the user did not set the channelBounds parameter.
            #set it empirically using the folloing hueristic.
            firstSampleForCB <- TRUE
            for (sampleName in samplesInExp) {
                #assumes that all exprsMat is in activeChannel column order.
                #this is determined by extractDataFromGS.
                exprsMat <- readRDS(file.path(normalizePath(projectPath),
                                              "faustData",
                                              "sampleData",
                                              sampleName,
                                              "exprsMat.rds"))
                if (firstSampleForCB) {
                    quantileDataLow <- apply(exprsMat,2,function(x){quantile(x,probs=0.01)})
                    quantileDataHigh <- apply(exprsMat,2,function(x){quantile(x,probs=0.99)})
                    empiricalCB <- matrix(-Inf,nrow=2,ncol=ncol(exprsMat))
                    colnames(empiricalCB) <- colnames(exprsMat)
                    rownames(empiricalCB) <- c("Low","High")
                    firstSampleForCB <- FALSE
                }
                else {
                    quantileDataLow <- rbind(quantileDataLow,
                                             apply(exprsMat,2,function(x){quantile(x,probs=0.01)}))
                    quantileDataHigh <- rbind(quantileDataHigh,
                                              apply(exprsMat,2,function(x){quantile(x,probs=0.99)}))
                }
            }
            empiricalCB["Low",] <- as.numeric(apply(quantileDataLow,2,function(x){quantile(x,probs=0.05)}))
            empiricalCB["High",] <- as.numeric(apply(quantileDataHigh,2,function(x){quantile(x,probs=0.95)}))
            print("Set the channelBounds parameter empirically. Resulted in the following matrix:")
            print(empiricalCB)
            if (uniqueHNum==1) {
                channelBoundsUsedByFAUST <- empiricalCB
            }
            else {
                internalCBList <- list()
                for (hName in uniqueHierarchyNames) {
                    internalCBList <- append(internalCBList,list(empiricalCB))
                    names(internalCBList)[length(internalCBList)] <- hName
                }
                channelBoundsUsedByFAUST <- internalCBList
            }
        }
        else if ((is.matrix(channelBounds)) && (uniqueHNum > 1)) {
            internalCBList <- list()
            for (hName in uniqueHierarchyNames) {
                internalCBList <- append(internalCBList,list(channelBounds))
                names(internalCBList)[length(internalCBList)] <- hName
            }
            channelBoundsUsedByFAUST <- internalCBList
        }
        saveRDS(channelBoundsUsedByFAUST,
                file.path(normalizePath(projectPath),
                          "faustData",
                          "metaData",
                          "channelBoundsUsedByFAUST.rds"))
        for (sampleName in samplesInExp) {
            exprsMat <- readRDS(file.path(normalizePath(projectPath),
                                          "faustData",
                                          "sampleData",
                                          sampleName,
                                          "exprsMat.rds"))
            resMat <- matrix(0,nrow = nrow(exprsMat), ncol = ncol(exprsMat))
            colnames(resMat) <- colnames(exprsMat)
            if (uniqueHNum > 1) {
                currentHID <- analysisMap[which(analysisMap$sampleName ==sampleName),"impH"]
                currentChannelBounds <- channelBoundsUsedByFAUST[[currentHID]]
            }
            else {
                #we assume that the channel bounds are a matrix in the event there is only
                #one imputation level.
                currentChannelBounds <- channelBoundsUsedByFAUST
            }
            if (debugFlag) {
                print("Extracting from sample: ")
                print(sampleName)
                print("Using restriction bounds: ")
                print(currentChannelBounds)
            }
            for (channel in colnames(currentChannelBounds)) {
                lowLookup <- which(exprsMat[,channel] <= currentChannelBounds["Low",channel])#nolint
                highLookup <- which(exprsMat[,channel] >= currentChannelBounds["High",channel])#nolint
                if (length(lowLookup)) {
                    resMat[lowLookup,channel] <- 1
                }
                if (length(highLookup)) {
                    resMat[highLookup,channel] <- 2
                }
            }
            saveRDS(resMat,
                    file.path(normalizePath(projectPath),
                              "faustData",
                              "sampleData",
                              sampleName,
                              "resMat.rds"))
        }
        madeResMats <- TRUE
        saveRDS(madeResMats,
                file.path(normalizePath(projectPath),
                          "faustData",
                          "metaData",
                          "madeResMats.rds"))
    }
    return()
}    
