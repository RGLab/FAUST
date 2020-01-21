.processChannelBounds <- function(
                                  samplesInExp,
                                  projectPath,
                                  channelBounds,
                                  analysisMap,
                                  debugFlag
                                  )
{
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
        if (debugFlag) {
            print("Set the channelBounds parameter empirically. Resulted in the following matrix:")
            print(empiricalCB)
        }
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
    #saveRDS(channelBoundsUsedByFAUST,
    #        file.path(normalizePath(projectPath),
    #                  "faustData",
    #                  "metaData",
    #                  "channelBoundsUsedByFAUST.rds"))
    #test to see if the channel bounds have changed.
    if (!file.exists(file.path(normalizePath(projectPath),
                               "faustData",
                               "metaData",
                               "channelBoundsUsedByFAUST.rds")))
    {
#                               "channelBounds.rds"))) {
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
            uniqueLevels <- unique(analysisMap[,"analysisLevel"])
            #delete flags for levels with annotation forests in order to regrow them.
            #also delete flags for scamp clsuterings to relabel data.
            for (analysisLevel in uniqueLevels) {
                if (file.exists(file.path(normalizePath(projectPath),
                                          "faustData",
                                          "levelData",
                                          analysisLevel,
                                          "aLevelComplete.rds"))) {
                    file.remove(file.path(normalizePath(projectPath),
                                          "faustData",
                                          "levelData",
                                          analysisLevel,
                                          "aLevelComplete.rds"))
                }
                if (file.exists(file.path(normalizePath(projectPath),
                                          "faustData",
                                          "levelData",
                                          analysisLevel,
                                          "scampALevelComplete.rds"))) {
                    file.remove(file.path(normalizePath(projectPath),
                                          "faustData",
                                          "levelData",
                                          analysisLevel,
                                          "scampALevelComplete.rds"))
                }
            }
            #if (file.exists(file.path(normalizePath(projectPath),
            #                          "faustData",
            #                          "metaData",
            #                          "parsedGS.rds"))) {
            #    file.remove(file.path(normalizePath(projectPath),
            #                          "faustData",
            #                          "metaData",
            #                          "parsedGS.rds"))
            #}
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
    
