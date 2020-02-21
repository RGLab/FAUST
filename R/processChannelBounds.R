.processChannelBounds <- function(samplesInExp,
                                  projectPath,
                                  channelBounds,
                                  debugFlag)
{
    analysisMap <- readRDS(file.path(normalizePath(projectPath),
                                     "faustData",
                                     "metaData",
                                     "analysisMap.rds"))

    uniqueHierarchyNames <- names(table(analysisMap[,"impH",drop=TRUE]))
    uniqueHNum <- length(uniqueHierarchyNames)
    channelBoundsUsedByFAUST <- channelBounds
    if ((is.character(channelBounds)) && (channelBounds=="")) {
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
                #initialize containers for channel bounds data

                #first percentile matrix
                quantileDataLowPrep <- apply(exprsMat,2,
                                             function(x){quantile(x,probs=0.01)})
                quantileDataLow <- matrix(quantileDataLowPrep,
                                          ncol=length(quantileDataLowPrep))
                colnames(quantileDataLow) <- names(quantileDataLowPrep)

                #99th percentile matrix
                quantileDataHighPrep <- apply(exprsMat,2,function(x){quantile(x,probs=0.99)})
                quantileDataHigh <- matrix(quantileDataHighPrep,
                                           ncol=length(quantileDataHighPrep))
                colnames(quantileDataHigh) <- names(quantileDataHighPrep)

                #container for empirical channel bounds
                empiricalCB <- matrix(-Inf,nrow=2,ncol=ncol(exprsMat))
                colnames(empiricalCB) <- colnames(exprsMat)
                rownames(empiricalCB) <- c("Low","High")
                firstSampleForCB <- FALSE
            }
            else {
                #update percentile matrices with the new dataset
                quantileDataLow <- rbind(quantileDataLow,
                                         apply(exprsMat,2,function(x){quantile(x,probs=0.01)}))
                quantileDataHigh <- rbind(quantileDataHigh,
                                          apply(exprsMat,2,function(x){quantile(x,probs=0.99)}))
            }
        }
        empiricalCB["Low",] <- as.numeric(apply(quantileDataLow,2,function(x){quantile(x,probs=0.05)}))
        empiricalCB["High",] <- as.numeric(apply(quantileDataHigh,2,function(x){quantile(x,probs=0.95)}))
        if (debugFlag) {
            print("Set the channelBounds parameter empirically. Resulted in the following matrix:")
            print(empiricalCB)
        }
        #if (uniqueHNum==1) {
        #    channelBoundsUsedByFAUST <- empiricalCB
        #}
        #else {

        #switch to always using a list.
        internalCBList <- list()
        for (hName in uniqueHierarchyNames) {
            internalCBList <- append(internalCBList,list(empiricalCB))
            names(internalCBList)[length(internalCBList)] <- hName
        }
        channelBoundsUsedByFAUST <- internalCBList
    }
    else if (is.matrix(channelBounds)) { #&& (uniqueHNum > 1)) {
        internalCBList <- list()
        for (hName in uniqueHierarchyNames) {
            internalCBList <- append(internalCBList,list(channelBounds))
            names(internalCBList)[length(internalCBList)] <- hName
        }
        channelBoundsUsedByFAUST <- internalCBList
    }
    #test to see if the channel bounds have changed.
    if (!file.exists(file.path(normalizePath(projectPath),
                               "faustData",
                               "metaData",
                               "channelBoundsUsedByFAUST.rds")))
    {
        #write the channel bounds provided at the faust call to disk.
        saveRDS(channelBoundsUsedByFAUST,file.path(normalizePath(projectPath),
                                                   "faustData",
                                                   "metaData",
                                                   "channelBoundsUsedByFAUST.rds"))
    }
    else {
        #if the channel bounds exist already, the faust method has been run at least once.
        #check to see if any modifications have been made to the channel bounds.
        oldChannelBounds <- readRDS(file.path(normalizePath(projectPath),
                                              "faustData",
                                              "metaData",
                                              "channelBoundsUsedByFAUST.rds"))
        if ((!identical(oldChannelBounds,channelBoundsUsedByFAUST)) &&
            (file.exists(file.path(normalizePath(projectPath),
                                   "faustData",
                                   "metaData",
                                   "bigForestDone.rds")))) {
            if (debugFlag) print("Detected change to channelBoundsUsedByFAUST.")
            file.remove(file.path(normalizePath(projectPath),
                                  "faustData",
                                  "metaData",
                                  "bigForestDone.rds"))
            if (file.exists(file.path(normalizePath(projectPath),
                                      "faustData",
                                      "metaData",
                                      "channelBounds.rds"))) {
                file.remove(file.path(normalizePath(projectPath),
                                      "faustData",
                                      "metaData",
                                      "channelBounds.rds"))
                saveRDS(channelBounds,
                        file.path(normalizePath(projectPath),
                                  "faustData",
                                  "metaData",
                                  "channelBounds.rds"))
            }
            file.remove(file.path(normalizePath(projectPath),
                                  "faustData",
                                  "metaData",
                                  "channelBoundsUsedByFAUST.rds"))
            file.remove(file.path(normalizePath(projectPath),
                                  "faustData",
                                  "metaData",
                                  "madeResMats.rds"))
            uniqueExpUnits <- unique(analysisMap[,"experimentalUnit"])
            #delete bools for experimental units that already have
            #annotation forests in order to regrow them.
            #also delete bools for scamp clusterings to relabel data.
            for (experimentalUnit in uniqueExpUnits) {
                if (file.exists(file.path(normalizePath(projectPath),
                                          "faustData",
                                          "expUnitData",
                                          experimentalUnit,
                                          "expUnitComplete.rds"))) {
                    file.remove(file.path(normalizePath(projectPath),
                                          "faustData",
                                          "expUnitData",
                                          experimentalUnit,
                                          "expUnitComplete.rds"))
                }
                if (file.exists(file.path(normalizePath(projectPath),
                                          "faustData",
                                          "expUnitData",
                                          experimentalUnit,
                                          "scampALevelComplete.rds"))) {
                    file.remove(file.path(normalizePath(projectPath),
                                          "faustData",
                                          "expUnitData",
                                          experimentalUnit,
                                          "scampALevelComplete.rds"))
                }
            }
            if (file.exists(file.path(normalizePath(projectPath),
                                      "faustData",
                                      "metaData",
                                      "firstALReady.rds"))) {
                file.remove(file.path(normalizePath(projectPath),
                                      "faustData",
                                      "metaData",
                                      "firstALReady.rds"))
            }
            saveRDS(channelBoundsUsedByFAUST,file.path(normalizePath(projectPath),
                                                       "faustData",
                                                       "metaData",
                                                       "channelBoundsUsedByFAUST.rds"))
        }
    }
    return()
}

