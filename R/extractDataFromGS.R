.extractDataFromGS <- function(gs,
                               analysisMap,
                               activeChannels,
                               channelBounds,
                               startingCellPop,
                               projectPath=".") {
    if (!file.exists(paste0(projectPath,"/faustData/metaData/parsedGS.rds"))) {
        samplesInGS <- flowWorkspace::sampleNames(gs)
        for (sampleName in samplesInGS) {
            dataSet <- flowWorkspace::getData(gs[[sampleName]],startingCellPop)
            exprsMatIn <- flowCore::exprs(dataSet)
            markers <- Biobase::pData(flowCore::parameters(dataSet))
            colnames(exprsMatIn) <- as.vector(markers[match(colnames(exprsMatIn),markers[,"name"], nomatch=0),]$desc)#nolint
            exprsMat <- exprsMatIn[,activeChannels]
            if (!dir.exists(paste0(projectPath,"/faustData/sampleData/",sampleName))) {#nolint
                dir.create(paste0(projectPath,"/faustData/sampleData/",sampleName))
            }
            saveRDS(exprsMat,paste0(projectPath,"/faustData/sampleData/",sampleName,"/exprsMat.rds"))#nolint
            resMat <- matrix(0,nrow = nrow(exprsMat), ncol = ncol(exprsMat))
            colnames(resMat) <- colnames(exprsMat)
            for (channel in colnames(channelBounds)) {
                lowLookup <- which(exprsMat[,channel] <= channelBounds["Low",channel])#nolint
                highLookup <- which(exprsMat[,channel] >= channelBounds["High",channel])#nolint
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
