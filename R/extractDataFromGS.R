.extractDataFromGS <- function(gs,
                               activeChannels,
                               startingCellPop,
                               projectPath=".",
                               debugFlag
                               ) {
    if (!file.exists(paste0(projectPath,"/faustData/metaData/parsedGS.rds"))) {
        samplesInGS <- flowWorkspace::sampleNames(gs)
        for (sampleName in samplesInGS) {
            if (debugFlag) {
                print(paste0("Extracting data from sample: ",sampleName))
            }
            dataSet <- flowWorkspace::getData(gs[[sampleName]],startingCellPop)
            exprsMatIn <- flowCore::exprs(dataSet)
            markers <- Biobase::pData(flowCore::parameters(dataSet))
            colnames(exprsMatIn) <- as.vector(markers[match(colnames(exprsMatIn),markers[,"name"], nomatch=0),]$desc)
            exprsMat <- exprsMatIn[,activeChannels]
            if (!dir.exists(paste0(projectPath,"/faustData/sampleData/",sampleName))) {#nolint
                dir.create(paste0(projectPath,"/faustData/sampleData/",sampleName))
            }
            saveRDS(exprsMat,paste0(projectPath,"/faustData/sampleData/",sampleName,"/exprsMat.rds"))#nolint
        }
        parsedGS <- TRUE
        saveRDS(parsedGS,paste0(projectPath,"/faustData/metaData/parsedGS.rds"))
    }
    return()
}    
