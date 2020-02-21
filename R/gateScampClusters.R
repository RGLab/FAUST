.gateScampClusters <- function(projectPath,
                               debugFlag)
{


    analysisMap <- readRDS(file.path(normalizePath(projectPath),
                                     "faustData",
                                     "metaData",
                                     "analysisMap.rds"))

    startingCellPop <- readRDS(file.path(normalizePath(projectPath),
                                         "faustData",
                                         "metaData",
                                         "sanitizedCellPopStr.rds"))

    selectedChannels <- readRDS(file.path(normalizePath(projectPath),
                                          "faustData",
                                          "gateData",
                                          paste0(startingCellPop,"_selectedChannels.rds")))

    resList <- readRDS(file.path(normalizePath(projectPath),
                                 "faustData",
                                 "gateData",
                                 paste0(startingCellPop,"_resList.rds")))
    activeSamples <- analysisMap[,"sampleName"]
    scampCellPops <- readRDS(file.path(normalizePath(projectPath),
                                       "faustData",
                                       "metaData",
                                       "scampClusterNames.rds"))
    for (sampleName in activeSamples) {
        sAnn <- utils::read.table(file = file.path(normalizePath(projectPath),
                                                   "faustData",
                                                   "sampleData",
                                                   sampleName,
                                                   "scampAnnotation.csv"),
                                  header = FALSE, sep = "`",
                                  stringsAsFactors = FALSE)[,1]
        annotationMatrix <- utils::read.table(file = file.path(normalizePath(projectPath),
                                                               "faustData",
                                                               "sampleData",
                                                               sampleName,
                                                               "annotationMatrix.csv"),
                                       header = FALSE,
                                       sep = ",",
                                       stringsAsFactors = FALSE)
        colnames(annotationMatrix) <- selectedChannels
        if (!(length(sAnn) == nrow(annotationMatrix))) {
            print(paste0("Annotation matrix doesn't match scamp annotation in ", sampleName))
            stop("Investigate.")
        }
        expUnit <- analysisMap[which(analysisMap[,"sampleName"] == sampleName), "experimentalUnit"]
        scampAF <- rep(list(NA),length(selectedChannels))
        names(scampAF) <- selectedChannels
        for (channel in selectedChannels) {
            scampAF[[channel]] <- resList[[channel]][[expUnit]]
        }
        gateNums <- unlist(lapply(scampAF,length)) + 1
        exactPartition <- gateSample(as.matrix(annotationMatrix), selectedChannels, gateNums, scampCellPops);

        # exactPartition <- rep("0_0_0_0_0",nrow(annotationMatrix))
        # for (rowNum in seq(nrow(annotationMatrix))) {
        #     ep <- paste0(paste0(paste0(selectedChannels,"~",annotationMatrix[rowNum,]),"~",gateNums,"~"),collapse="")
        #     if (ep %in% scampCellPops) {
        #         exactPartition[rowNum] <- ep
        #     }
        # }
        data.table::fwrite(list(exactPartition),
                           file = file.path(normalizePath(projectPath),
                                            "faustData",
                                            "sampleData",
                                            sampleName,
                                            "faustAnnotation.csv"),
                           sep = "~",
                           append = FALSE,
                           row.names = FALSE,
                           col.names = FALSE,
                           quote = FALSE)
    }
    return()
}
