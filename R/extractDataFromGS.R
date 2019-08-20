.extractDataFromGS <- function(gs,
                               analysisMap,
                               activeChannels,
                               channelBounds,
                               startingCellPop,
                               projectPath=".",
                               debugFlag
                               ) {
    if (!file.exists(paste0(projectPath,"/faustData/metaData/parsedGS.rds"))) {
        samplesInGS <- flowWorkspace::sampleNames(gs)
        uniqueHNum <- length(table(analysisMap[,"impH",drop=TRUE]))
        for (sampleName in samplesInGS) {
            dataSet <- flowWorkspace::getData(gs[[sampleName]],startingCellPop)
            exprsMatIn <- flowCore::exprs(dataSet)
            markers <- Biobase::pData(flowCore::parameters(dataSet))
            colnames(exprsMatIn) <- as.vector(markers[match(colnames(exprsMatIn),markers[,"name"], nomatch=0),]$desc)
            exprsMat <- exprsMatIn[,activeChannels]
            if (!dir.exists(paste0(projectPath,"/faustData/sampleData/",sampleName))) {#nolint
                dir.create(paste0(projectPath,"/faustData/sampleData/",sampleName))
            }
            saveRDS(exprsMat,paste0(projectPath,"/faustData/sampleData/",sampleName,"/exprsMat.rds"))#nolint
            resMat <- matrix(0,nrow = nrow(exprsMat), ncol = ncol(exprsMat))
            colnames(resMat) <- colnames(exprsMat)
            if (uniqueHNum > 1) {
                currentHID <- analysisMap[which(analysisMap$sampleName ==sampleName),"impH"]
                currentChannelBounds <- channelBounds[[currentHID]]
            }
            else {
                #we assume that the channel bounds are a matrix in the event there is only
                #one imputation level.
                currentChannelBounds <- channelBounds
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
            saveRDS(resMat,paste0(projectPath,"/faustData/sampleData/",sampleName,"/resMat.rds"))
        }
        parsedGS <- TRUE
        saveRDS(parsedGS,paste0(projectPath,"/faustData/metaData/parsedGS.rds"))
    }
    return()
}    
