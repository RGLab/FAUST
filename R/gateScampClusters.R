.gateScampClusters <- function(startingCellPop,
                               analysisMap,
                               selectedChannels,
                               projectPath=".",
                               debugFlag=FALSE)
{
    resList <- readRDS(paste0(projectPath,"/faustData/gateData/",startingCellPop,"_resList.rds"))
    activeSamples <- analysisMap[,"sampleName"]
    scampCellPops <- readRDS(paste0(projectPath,"/faustData/metaData/scampClusterNames.rds"))
    for (sampleName in activeSamples) {
        sAnn <- utils::read.table(file = paste0(projectPath,"/faustData/sampleData/",sampleName,"/scampAnnotation.csv"),
                                  header = FALSE, sep = "`", 
                                  stringsAsFactors = FALSE)[,1]
        annotationMatrix <- utils::read.table(file = paste0(projectPath,"/faustData/sampleData/",sampleName,"/annotationMatrix.csv"),
                                       header = FALSE,
                                       sep = ",",
                                       stringsAsFactors = FALSE)
        colnames(annotationMatrix) <- selectedChannels
        if (!(length(sAnn) == nrow(annotationMatrix))) {
            print(paste0("Annotation matrix doesn't match scamp annotation in ", sampleName))
            stop("Investigate.")
        }
        aLevel <- analysisMap[which(analysisMap[,"sampleName"] == sampleName), "analysisLevel"]
        scampAF <- rep(list(NA),length(selectedChannels))
        names(scampAF) <- selectedChannels
        for (channel in selectedChannels) {
            scampAF[[channel]] <- resList[[channel]][[aLevel]]
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
                    file = paste0(projectPath, "/faustData/sampleData/", sampleName, "/faustAnnotation.csv"),
                    sep = "~",
                    append = FALSE,
                    row.names = FALSE,
                    col.names = FALSE,
                    quote = FALSE)
    }
    return()
}
