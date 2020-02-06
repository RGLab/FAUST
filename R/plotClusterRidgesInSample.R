#' Fast Annotation using Shape-constrained Trees
#'
#' This function plots the
#' strategy for a population found by the FAUST pipeline.
#'
#' @param clusterPhenotype A text string corresponding to a cell population
#' determined by FAUST. Must be taken from the
#' projectPath/faustData/faustCountMatrix.rds column names.
#'
#' @param sampleName A text string corresponding to the experimental
#' sample in which you want to plot clusterPhenotype ridgelines.
#' Must be a directory name in projectPath/faustData/sampleData.
#'
#' @param outputDirName A text string which will be create a directory
#' to store plot output. The directory will be placed in
#' projectPath/faustData/plotData/clusterRidges.
#'
#' @param projectPath An absolute path to the directory where you
#' ran the FAUST.
#'
#' @param plottingDevice string with device for saving graphical output.
#' By default it is set to "pdf".
#'
#' @return plotClusterRidgesInSample returns a null value on completion.
#' The main output is a file of type plottingDevice located in
#' "projectPath/faustData/plotData/clusterRidges/outputDirName"
#' that shows ridgelines of the clusterPhenotype distribution set against the
#' entire distribution of the sampleName for each marker used to define
#' the phenotype.
#'
#' @examples
#' \dontrun{
#' myProjectPath <- "/path/to/directory/housing/faustData"
#' countMatrix <- readRDS(paste0(myProjectPath,"/faustData/faustCountMatrix.rds"))
#' allSnames <- list.files(file.path(myProjectPath,
#'                                   "faustData",
#'                                   "sampleData"))
#'
#' plotClusterRidgesInSample(
#'     clusterPhenotype=colnames(countMatrix)[1],
#'     sampleName=allSnames[1],
#'     outputDirName=testPlot,
#'     projectPath=myProjectPath,
#'     plottingDevice="pdf"
#' )
#' }
#' @export
#' @md
plotClusterRidgesInSample <- function(
                                      clusterPhenotype,
                                      sampleName,
                                      outputDirName,
                                      projectPath,
                                      plottingDevice="pdf"
                                      )
{
    if (!dir.exists(file.path(normalizePath(projectPath),
                              "faustData"))) {
        print(paste0("Do not detect a faustData directory at the following projectPath."))
        print(projectPath)
        stop("Set projectPath to point to a completed faust run.")
    }
    if (!dir.exists(file.path(normalizePath(projectPath),
                              "faustData",
                              "plotData",
                              "clusterRidges")))
    {
        dir.create(file.path(normalizePath(projectPath),
                             "faustData",
                             "plotData",
                             "clusterRidges"))
    }
    cnMap <- readRDS(file.path(normalizePath(projectPath),
                               "faustData",
                               "metaData",
                               "colNameMap.rds"))
    clusterLookup <- which(cnMap[,"newColNames"]==clusterPhenotype)
    if (length(clusterLookup)==0) {
        print(paste0("The clusterPhenotype parameter is set to the following string."))
        print(clusterPhenotype)
        print("This cluster is not defined in the completed faust run.")
        stop("Set clusterPhenotype to a column value in faustCountMatrix.rds and run again.")
    }
    encodedCluster <- cnMap[clusterLookup,"faustColNames"]
    if (!dir.exists(file.path(normalizePath(projectPath),
                              "faustData",
                              "plotData",
                              "clusterRidges",
                              outputDirName)))
    {
        dir.create(file.path(normalizePath(projectPath),
                             "faustData",
                             "plotData",
                             "clusterRidges",
                             outputDirName))
    }
    allSamples <- list.files(file.path(normalizePath(projectPath),
                                       "faustData",
                                       "sampleData"))
    snLookup <- which(allSamples==sampleName)
    if (length(snLookup)==0) {
        print(paste0("sampleName not detected in ",
                     file.path(normalizePath(projectPath),
                               "faustData",
                               "sampleData")))
        stop("Set sampleName to a sub-directory located there.")
    }
    selCols <- readRDS(file.path(normalizePath(projectPath),
                                 "faustData",
                                 "metaData",
                                 "initSelC.rds"))
    exprsMatIn <- readRDS(file.path(normalizePath(projectPath),
                                    "faustData",
                                    "sampleData",
                                    sampleName,
                                    "exprsMat.rds"))
    sampleExprs <- as.data.frame(exprsMatIn[,selCols,drop=FALSE])
    fa <- read.table(file=file.path(normalizePath(projectPath),
                                    "faustData",
                                    "sampleData",
                                    sampleName,
                                    "faustAnnotation.csv"),
                     header=FALSE,
                     sep="`")[,1]
    plotDFPrep <- gather(sampleExprs,key="Channel",value="Expression")
    plotDFPrep$Source <- "Sample\nDistribution"
    cluster1Lookup <- which(fa==encodedCluster)
    if (length(cluster1Lookup)) {
        clusterExprs <- as.data.frame(sampleExprs[cluster1Lookup,,drop=FALSE])
        clusterPlotDF <- gather(clusterExprs,key="Channel",value="Expression")
        clusterPlotDF$Source <- "FAUST Cluster\nDistribution"
        plotDFPrep2 <- rbind(plotDFPrep,clusterPlotDF)
        qs <- as.numeric(quantile(plotDFPrep2[,"Expression"],probs=c(0.01,0.99)))
        plotLookup <- intersect(which(plotDFPrep2[,"Expression"] >= qs[1]),
                                which(plotDFPrep2[,"Expression"] <= qs[2]))
        plotDF <- plotDFPrep2[plotLookup,,drop=FALSE]
        pRidges <- ggplot(plotDF,aes_string(x="Expression",y="Source",fill="Source"))+
            geom_density_ridges(scale=1,jittered_points = TRUE, point_size = 1.75, point_shape = "|",
                                position = position_points_jitter(width = 0.01, height = 0))+
            scale_fill_viridis(begin=0.25,end=0.75,discrete=TRUE)+
            facet_wrap(~Channel,nrow=2,scales="free")+
            theme_bw(base_size = 18)+
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1),
                legend.position="bottom",
                plot.title = element_text(hjust = 0.5),
                strip.text = element_text(size = 18)
            )+
            ggtitle("Sample and cluster expression distributions\nTruncated at 1st and 99th percentiles")+
            theme(axis.title.y=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank())
        ggsave(
            filename=file.path(normalizePath(projectPath),
                               "faustData",
                               "plotData",
                               "clusterRidges",
                               outputDirName,
                               paste0(sampleName,".",plottingDevice)),
            plot=pRidges,
            width=14,
            height=7,
            units="in",
            dpi=300
        )
    }
    else {
        print(paste0("Cluster is exactly zero in sample ",sampleName))
        print("No plot created.")
    }
    return()
}


