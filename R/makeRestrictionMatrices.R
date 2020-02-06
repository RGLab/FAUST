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
        channelBoundsUsedByFAUST <- readRDS(file.path(normalizePath(projectPath),
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

            #Extract the channel bounds matrix associated with the level of the
            #imputation hierarchy
            currentHID <- analysisMap[which(analysisMap$sampleName ==sampleName),"impH"]
            currentChannelBounds <- channelBoundsUsedByFAUST[[currentHID]]
            if (debugFlag) {
                print("Extracting from sample: ")
                print(sampleName)
                print("Using restriction bounds: ")
                print(currentChannelBounds)
            }
            #Construct the restriction matrix using the channel bounds
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
