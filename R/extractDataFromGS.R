.extractDataFromGS <- function(gs,
                               activeChannels,
                               startingCellPop,
                               projectPath,
                               debugFlag
                               )
{
    if (!file.exists(file.path(normalizePath(projectPath),
                               "faustData",
                               "metaData",
                               "parsedGS.rds"))) {
        samplesInGS <- flowWorkspace::sampleNames(gs)
        for (sampleName in samplesInGS) {
            if (debugFlag) {
                print(paste0("Extracting data from sample: ",sampleName))
            }
            #dataSet <- flowWorkspace::getData(gs[[sampleName]],startingCellPop)
            dataSet <- flowWorkspace::gh_pop_get_data(gs[[sampleName]],startingCellPop)
            exprsMatIn <- flowCore::exprs(dataSet)
            markers <- Biobase::pData(flowCore::parameters(dataSet))
            colnames(exprsMatIn) <- as.vector(markers[match(colnames(exprsMatIn),markers[,"name"], nomatch=0),]$desc)
            exprsMat <- exprsMatIn[,activeChannels,drop=FALSE]
            if (!dir.exists(file.path(normalizePath(projectPath),
                                      "faustData",
                                      "sampleData",
                                      sampleName))) {
                dir.create(file.path(normalizePath(projectPath),
                                     "faustData",
                                     "sampleData",
                                     sampleName))
            }
            saveRDS(exprsMat,
                    file.path(normalizePath(projectPath),
                              "faustData",
                              "sampleData",
                              sampleName,
                              "exprsMat.rds"))
        }
        parsedGS <- TRUE
        saveRDS(parsedGS,
                file.path(normalizePath(projectPath),
                          "faustData",
                          "metaData",
                          "parsedGS.rds"))
    }
    return()
}
