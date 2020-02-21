.mkAnnMats <- function(projectPath)
{
    parentNode <- readRDS(file.path(normalizePath(projectPath),
                                    "faustData",
                                    "metaData",
                                    "sanitizedCellPopStr.rds"))

    analysisMap <- readRDS(file.path(normalizePath(projectPath),
                                     "faustData",
                                     "metaData",
                                     "analysisMap.rds"))

    resList <- readRDS(file.path(normalizePath(projectPath),
                                 "faustData",
                                 "gateData",
                                 paste0(parentNode,"_resList.rds")))

    selC <- readRDS(file.path(normalizePath(projectPath),
                              "faustData",
                              "gateData",
                              paste0(parentNode,"_selectedChannels.rds")))
    for (sampleNum in seq(nrow(analysisMap))) {
        sampleName <- analysisMap[sampleNum,"sampleName"]
        expUnit <- analysisMap[sampleNum,"experimentalUnit"]
        exprsMatIn <- readRDS(file.path(normalizePath(projectPath),
                                        "faustData",
                                        "sampleData",
                                        sampleName,
                                        "exprsMat.rds"))
        exprsMat <- exprsMatIn[,selC,drop=FALSE]
        annotationMatrix <- matrix(0,nrow=nrow(exprsMat),ncol=ncol(exprsMat))
        colnames(annotationMatrix) <- colnames(exprsMat)
        for (column in names(resList)) {
            annotationMatrix[,column] <- 1
            gateVals <- resList[[column]][[expUnit]]
            if (is.list(gateVals)) {
                gateVals <- as.vector(unlist(gateVals))
            }
            for (gateVal in gateVals) {
                annLook <- which(exprsMat[,column] >= gateVal)
                if (length(annLook)) {
                    annotationMatrix[annLook,column] <- (annotationMatrix[annLook,column]+1)
                }
            }
        }
        data.table::fwrite(as.data.frame(annotationMatrix),
                           file=file.path(normalizePath(projectPath),
                                          "faustData",
                                          "sampleData",
                                          sampleName,
                                          "annotationMatrix.csv"),
                           sep=",",
                           row.names=FALSE,
                           col.names=FALSE)
    }
}

