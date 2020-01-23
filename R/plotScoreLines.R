.plotScoreLines <- function(
                            projectPath,
                            depthScoreThreshold,
                            selectionQuantile,
                            forceList,
                            plottingDevice
                            )
{
    scoreDF <- as.data.frame(readRDS(file.path(normalizePath(projectPath),
                                               "faustData",
                                               "metaData",
                                               "scoreMat.rds")))
    qScore <- as.data.frame(apply(scoreDF,2,
                                  function(x){as.numeric(stats::quantile(x,probs=seq(0,1,by=0.01)))}))
    qScore$Quantile <- seq(0,1,by=0.01)
    activeScore <- apply(scoreDF,2,function(x){as.numeric(stats::quantile(x,probs=selectionQuantile))})
    lineTypeMap <- rep("twodash",length(colnames(scoreDF)))
    names(lineTypeMap) <- colnames(scoreDF)
    solidNames <- names(which(activeScore >= depthScoreThreshold))
    if (length(solidNames)) lineTypeMap[solidNames] <- "solid"
    colorMap <- viridis::magma(length(colnames(scoreDF)))
    names(colorMap) <- colnames(scoreDF)
    selectedNames <- names(which(activeScore >= depthScoreThreshold))
    if (length(selectedNames)) colorMap[selectedNames] <- viridis::viridis(length(selectedNames))
    unselectedNames <- names(which(activeScore < depthScoreThreshold))
    if (length(unselectedNames)) colorMap[unselectedNames] <- viridis::viridis(length(unselectedNames))
    if (length(forceList)) {
        colorMap[names(forceList)] <- "#ff0000"
        lineTypeMap[names(forceList)] <- "solid"
    }
    qsDF <- tidyr::gather(qScore,key="Channel",value="QuantileValue",-Quantile)
    p <- ggplot(qsDF, aes(x = Quantile, y = QuantileValue, color = Channel)) +
      geom_line(aes(color = Channel, linetype = Channel)) +
      scale_color_manual(values = colorMap) +
      scale_linetype_manual(values = lineTypeMap) +
      geom_vline(
        xintercept = selectionQuantile,
        color = "red",
        linetype = "dotted",
        size = 0.5
      ) +
      geom_hline(
        yintercept = depthScoreThreshold,
        color = "red",
        linetype = "dotted",
        size = 0.5
      ) +
      xlab("Selection Quantile") +
      ylab("Channel Depth Score at Selection Quantile") +
        theme_bw()+
        theme(legend.position="bottom")
    fpNameOut <- file.path(normalizePath(projectPath),
                           "faustData",
                           "plotData",
                           paste0("scoreLines.",plottingDevice))
    ggplot2::ggsave(
                 filename=fpNameOut,
                 plot=p,
                 units="in",
                 height = 11,
                 width = 10
             )
    
    return()
}
