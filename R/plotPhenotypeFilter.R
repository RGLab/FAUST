.plotPhenotypeFilter <- function(
                                 projectPath,
                                 nameOccuranceNum,
                                 plottingDevice
                                 )
{
    
    nameSummary <- readRDS(file.path(normalizePath(projectPath),
                                     "faustData",
                                     "metaData",
                                     "scampNameSummary.rds"))
    nameSummaryPlotDF <- data.frame(x=seq(max(nameSummary)),
                                    y=sapply(seq(max(nameSummary)),function(x){
                                        length(which(nameSummary >= x))}))
    elbowLoc <- readRDS(
        file.path(normalizePath(projectPath),
                  "faustData",
                  "metaData",
                  "phenotypeElbowValue.rds")
    )
    nspOut <- ggplot(nameSummaryPlotDF,aes(x=x,y=y))+
        geom_line()+
        theme_bw()+
        geom_vline(xintercept=elbowLoc,col="red")+
        xlab("Number of times a phenotype appears across clusterings")+
        ylab("Number of phenotypes exceeding the appearance number")+
        ggtitle("Red line is nameOccuranceNum setting in faust")

    fpNameOut <- file.path(normalizePath(projectPath),
                           "faustData",
                           "plotData",
                           paste0("scampNamesPlot.",plottingDevice))
    ggplot2::ggsave(
                 filename=fpNameOut,
                 plot=nspOut,
                 units="in",
                 height = 8,
                 width = 8
             )
    return()
}


