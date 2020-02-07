.plotSampleHistograms <- function(
                                  projectPath,
                                  plottingDevice
                                  )
{
    analysisMap <- readRDS(file.path(normalizePath(projectPath),
                                     "faustData",
                                     "metaData",
                                     "analysisMap.rds"))

    startingCellPop <- readRDS(file.path(normalizePath(projectPath),
                                         "faustData",
                                         "metaData",
                                         "sanitizedCellPopStr.rds"))

    resList <- readRDS(file.path(normalizePath(projectPath),
                                 "faustData",
                                 "gateData",
                                 paste0(startingCellPop,"_resList.rds")))
    selC <- readRDS(file.path(normalizePath(projectPath),
                              "faustData",
                              "gateData",
                              paste0(startingCellPop,"_selectedChannels.rds")))

    for (sampleName in analysisMap[,"sampleName"]) {
        exprsMat <- readRDS(file.path(normalizePath(projectPath),
                                      "faustData",
                                      "sampleData",
                                      sampleName,
                                      "exprsMat.rds"))
        aLevel <- analysisMap[which(analysisMap[,"sampleName"]==sampleName),"analysisLevel"]
        plotList <- list()
        #initialize directories for the selected channels
        for (channel in selC) {
            sanitizedChannel <- gsub("[[:punct:]]","",channel)
            sanitizedChannel <- gsub("[[:space:]]","",sanitizedChannel)
            sanitizedChannel <- gsub("[[:cntrl:]]","",sanitizedChannel)
            if (!dir.exists(file.path(normalizePath(projectPath),
                                      "faustData",
                                      "plotData",
                                      "histograms",
                                      sanitizedChannel)))
            {
                dir.create(file.path(normalizePath(projectPath),
                                     "faustData",
                                     "plotData",
                                     "histograms",
                                     sanitizedChannel))
            }
        }
        #plot marker histograms by sample
        for (channel in selC) {
            channelData <- as.data.frame(exprsMat[,channel,drop=FALSE])
            colnames(channelData) <- "x"
            gateData <- resList[[channel]][[aLevel]]
            channelQs <- as.numeric(quantile(channelData$x,probs=c(0.01,0.99)))
            histLookupLow <- which(channelData$x >= channelQs[1])
            histLookupHigh <- which(channelData$x <= channelQs[2])
            histLookup <- intersect(histLookupLow,histLookupHigh)
            histData <- channelData[histLookup,"x",drop=FALSE]
            p <- .getHistogram(histData,channel,gateData)
            sanitizedChannel <- gsub("[[:punct:]]","",channel)
            sanitizedChannel <- gsub("[[:space:]]","",sanitizedChannel)
            sanitizedChannel <- gsub("[[:cntrl:]]","",sanitizedChannel)
            fpNameOut <- file.path(normalizePath(projectPath),
                                   "faustData",
                                   "plotData",
                                   "histograms",
                                   sanitizedChannel,
                                   paste0(sampleName,".",plottingDevice))
            ggplot2::ggsave(
                         filename=fpNameOut,
                         plot=p,
                         units="in",
                         height = 6,
                         width = 6
                     )
        }
    }
    return()
}

.getHistogram <- function(histData,channelName,gates) {
    fdBreaks <- pretty(range(histData[,"x"]),
                     n = grDevices::nclass.FD(histData[,"x"]), min.n = 1)
    binWidth <- fdBreaks[2]-fdBreaks[1]
    p <- ggplot(histData,aes(x=x)) +
        geom_histogram(binwidth=binWidth) +
        theme_bw()+
        geom_vline(xintercept=gates,color="red",linetype="dashed")+
        ggtitle(channelName)
    return(p)
}

