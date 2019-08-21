.plotMarkerHistograms <- function(analysisMap,
                                  startingCellPop,
                                  projectPath=".")
{
    resList <- readRDS(paste0(projectPath,"/faustData/gateData/",
                              startingCellPop,"_resList.rds"))
    selC <- readRDS(paste0(projectPath,"/faustData/gateData/",
                           startingCellPop,"_selectedChannels.rds"))
    uniqueLevels <- sort(unique(analysisMap[,"analysisLevel"]))
    firstALevel <- TRUE
    for (aLevel in uniqueLevels) {
        levelExprs <- readRDS(paste0(projectPath,"/faustData/levelData/",aLevel,"/levelExprs.rds"))
        rndLook <- sample(seq(nrow(levelExprs)),min(1000,nrow(levelExprs)))
        if (firstALevel) {
            plotMat <- levelExprs[rndLook,selC,drop=FALSE]
            firstALevel <- FALSE
        }
        else {
            plotMat <- rbind(plotMat,levelExprs[rndLook,selC,drop=FALSE])
        }
    }
    for (channel in selC) {
        channelData <- as.data.frame(plotMat[,channel,drop=FALSE])
        colnames(channelData) <- "x"
        channelQs <- as.numeric(quantile(channelData$x,probs=c(0.01,0.99)))
        histLookupLow <- which(channelData$x >= channelQs[1])
        histLookupHigh <- which(channelData$x <= channelQs[2])
        histLookup <- intersect(histLookupLow,histLookupHigh)
        histData <- channelData[histLookup,"x",drop=FALSE]
        gateData <- resList[[channel]]
        totalGateNum <- length(gateData[[1]])
        for (gateNum in seq(totalGateNum)) {
            gatesForPlotting <- unlist(lapply(gateData,function(x){x[gateNum]}))
            p <- .getHistogram(histData,channel,gatesForPlotting)
            cowplot::save_plot(paste0(projectPath,"/faustData/plotData/hist_",channel,"_ab_",gateNum,".pdf"),
                               p,
                               base_height = 7,
                               base_width = 7)
        }
    }
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
