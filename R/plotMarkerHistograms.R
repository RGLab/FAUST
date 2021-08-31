.plotMarkerHistograms <- function(
                                  projectPath,
                                  plottingDevice
                                  )
{
    #make aggregate histograms for diagnostics.
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

    uniqueExpUnits <- sort(unique(analysisMap[,"experimentalUnit"]))
    firstExpUnit <- TRUE
    for (expUnit in uniqueExpUnits) {
        expUnitExprs <- readRDS(file.path(normalizePath(projectPath),
                                        "faustData",
                                        "expUnitData",
                                        expUnit,
                                        "expUnitExprs.rds"))
        rndLook <- sample(seq(nrow(expUnitExprs)),min(1000,nrow(expUnitExprs)))
        if (firstExpUnit) {
            plotMat <- expUnitExprs[rndLook,selC,drop=FALSE]
            firstExpUnit <- FALSE
        }
        else {
            plotMat <- rbind(plotMat,expUnitExprs[rndLook,selC,drop=FALSE])
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
            fpNameOut <- file.path(normalizePath(projectPath),
                                   "faustData",
                                   "plotData",
                                   paste0("hist_",channel,"_ab_",gateNum,".",plottingDevice))

            ggplot2::ggsave(
                         filename=fpNameOut,
                         plot=p,
                         units="in",
                         height = 6,
                         width = 6
                     )
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
