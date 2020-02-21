#' Fast Annotation using Shape-constrained Trees
#'
#' This function plots the depth-score-order gating
#' strategy for a population found by the FAUST pipeline.
#'
#' @param cellPop A text string corresponding to a cell population
#' determined by the FAUST pipeline. Must be taken from the
#' projectPath/faustData/faustCountMatrix.rds column names.
#'
#' @param sampleName A text string corresponding to the experimental
#' sample in which you want to plot cellPop's gating strategy. Must be
#' taken from the projectPath/faustData/faustCountMatrix.rds row names.
#'
#' @param projectPath An absolute path to the directory where you
#' ran the FAUST pipeline.
#'
#' @return plotFaustGates returns a null value on completion. The main output is a pdf
#' located in
#' "projectPath/faustData/plotData/gatingStrats/cellPops"
#' that shows the gating strategy.
#'
#' @examples
#' \dontrun{
#' myProjectPath <- "/path/to/directory/housing/faustData"
#' countMatrix <- readRDS(paste0(myProjectPath,"/faustData/faustCountMatrix.rds"))
#' plotFaustGates(rownames(countMatrix)[1],colnames(countMatrix)[1],myProjectPath)
#' }
#' @export
#' @md
plotFaustGates <- function(cellPop,sampleName,projectPath=".") {
    if (cellPop=="0_0_0_0_0") return(NA)
    if (!dir.exists(file.path(normalizePath(projectPath),
                              "faustData")))
    {
        print(paste0("faustData directory not detected in ",projectPath))
        print("Update the projectPath variable.")
        return(NA)
    }
    if (!dir.exists(file.path(normalizePath(projectPath),
                              "faustData",
                              "plotData")))
    {
        print(paste0("faustData/plotData directory not detected in ",projectPath))
        print("This indicates the FAUST pipeline has not been run.")
        print("If it has, update the projectPath to the directory containing faustData")
        print("If it has not, run the FAUST pipeline before calling this function.")
        return(NA)
    }
    if (!dir.exists(file.path(normalizePath(projectPath),
                              "faustData",
                              "plotData",
                              "gatingStrats")))
    {
        dir.create(file.path(normalizePath(projectPath),
                             "faustData",
                             "plotData",
                             "gatingStrats"))
    }
    if (!dir.exists(file.path(normalizePath(projectPath),
                              "faustData",
                              "plotData",
                              "gatingStrats",
                              cellPop)))
    {
        dir.create(file.path(normalizePath(projectPath),
                             "faustData",
                             "plotData",
                             "gatingStrats",
                             cellPop))
    }
    colNameMap <- readRDS(file.path(normalizePath(projectPath),
                                    "faustData",
                                    "metaData",
                                    "colNameMap.rds"))
    cellPopLookup <- which(colNameMap[,"newColNames"]==cellPop)
    if (length(cellPopLookup)) {
        faustPopName <- colNameMap[cellPopLookup,"faustColNames"]
        p <- .plotFaustGates(cellPop=faustPopName,sampleName=sampleName,projectPath=projectPath)
        cowplot::save_plot(file.path(normalizePath(projectPath),
                                     "faustData",
                                     "plotData",
                                     "gatingStrats",
                                     cellPop,
                                     paste0(sampleName,".pdf")),
                           p,
                           base_height = 20,
                           base_width = 20)
        return()
    }
    else {
        return(NA)
    }
}

.plotFaustGates <- function(cellPop,sampleName,projectPath) {
    channelBounds <- readRDS(file.path(normalizePath(projectPath),
                                       "faustData",
                                       "metaData",
                                       "channelBounds.rds"))
    analysisMap <- readRDS(file.path(normalizePath(projectPath),
                                     "faustData",
                                     "metaData",
                                     "analysisMap.rds"))
    startingCellPop <- readRDS(file.path(normalizePath(projectPath),
                                         "faustData",
                                         "metaData",
                                         "startingCellPop.rds"))
    startingCellPop <- gsub("[[:punct:]]","",startingCellPop)
    startingCellPop <- gsub("[[:space:]]","",startingCellPop)
    startingCellPop <- gsub("[[:cntrl:]]","",startingCellPop)
    resList <- readRDS(file.path(normalizePath(projectPath),
                                 "faustData",
                                 "gateData",
                                 paste0(startingCellPop,"_resList.rds")))
    staticPop <- cellPop
    snLookup <- which(analysisMap[,"sampleName"] == sampleName)
    if (length(snLookup)) {
        expUnit <- analysisMap[snLookup,"experimentalUnit"]
        exprsMat <- readRDS(file.path(normalizePath(projectPath),
                                      "faustData",
                                      "sampleData",
                                      sampleName,
                                      "exprsMat.rds"))
        indexRows <- rep(TRUE,nrow(exprsMat))
        targetPop <- staticPop
        plotPop <- strsplit(targetPop,"~")[[1]]
        plotList <- list()
        while(length(plotPop)) {
            nc1 <- plotPop[1]
            ex1 <- plotPop[2]
            ux1 <- plotPop[3]
            plotPop <- plotPop[-c(1,2,3)]
            if (length(plotPop)) {
                nc2 <- plotPop[1]
                ex2 <- plotPop[2]
                ux2 <- plotPop[3]
                plotPop <- plotPop[-c(1,2,3)]
                activeRows <- which(indexRows == TRUE)
                plotData <- as.data.frame(exprsMat[activeRows,c(nc1,nc2),drop=FALSE])
                c1qs <- as.numeric(quantile(exprsMat[,nc1],probs=c(0.01,0.99)))
                c2qs <- as.numeric(quantile(exprsMat[,nc2],probs=c(0.01,0.99)))
                gv1 <- resList[[nc1]][[expUnit]]
                gv2 <- resList[[nc2]][[expUnit]]
                aMat <- matrix(1,nrow=nrow(plotData),ncol=ncol(plotData))
                colnames(aMat) <- colnames(plotData)
                for (gv in gv1) {
                    annLook <- which(plotData[,nc1] >= gv)
                    if (length(annLook)) {
                        aMat[annLook,nc1] <- (aMat[annLook,nc1] + 1)
                    }
                }
                for (gv in gv2) {
                    annLook <- which(plotData[,nc2] >= gv)
                    if (length(annLook)) {
                        aMat[annLook,nc2] <- (aMat[annLook,nc2] + 1)
                    }
                }
                look1 <- which(aMat[,nc1] == as.numeric(ex1))
                look2 <- which(aMat[,nc2] == as.numeric(ex2))
                lookup <- sort(intersect(look1,look2))
                if (nrow(plotData)) {
                    colnames(plotData) <- c("x","y")
                    if (ex1 == "1") pex1 <- "- "
                    else if (ex1 == ux1) pex1 <- "+ "
                    else pex1 <- paste0("_",ex1,"_",ux1," ")
                    if (ex2 == "1") pex2 <- "- "
                    else if (ex2 == ux2) pex2 <- "+ "
                    else pex2 <- paste0("_",ex2,"_",ux2," ")
                    #drop outliers
                    nc1qs <- as.numeric(quantile(plotData[,"x"],probs=c(0.01,0.99)))
                    nc2qs <- as.numeric(quantile(plotData[,"y"],probs=c(0.01,0.99)))
                    di1 <- which(plotData[,"x"] <= nc1qs[1])
                    di2 <- which(plotData[,"x"] >= nc1qs[2])
                    di3 <- which(plotData[,"y"] <= nc2qs[1])
                    di4 <- which(plotData[,"y"] >= nc2qs[2])
                    dropIndices <- sort(unique(c(di1,di2,di3,di4)))
                    plotData <- plotData[-dropIndices,,drop=FALSE]
                    p <- ggplot(plotData,aes(x=x,y=y)) + geom_hex(bins=128) + theme_bw(base_size=16)+
                        xlab(nc1)+ylab(nc2)+
                        geom_vline(xintercept=gv1,linetype="dashed",color="red") +
                        geom_hline(yintercept=gv2,linetype="dashed",color="red") +
                        theme(legend.position="none")+
                        coord_cartesian(xlim=c(c1qs[1],c1qs[2]),ylim=c(c2qs[1],c2qs[2]))+
                        ggtitle(paste0(nc1,pex1," ",nc2,pex2,": ",
                                       length(lookup),"/",nrow(plotData)))
                    plotList <- append(plotList,list(p))
                }
                indexRows <- rep(FALSE,nrow(exprsMat))
                updateRows <- activeRows[lookup]
                indexRows[updateRows] <- TRUE
            }
            else {
                activeRows <- which(indexRows == TRUE)
                plotData <- as.data.frame(exprsMat[activeRows,c(nc1),drop=FALSE])
                c1qs <- as.numeric(quantile(exprsMat[,nc1],probs=c(0.01,0.99)))
                gv1 <- resList[[nc1]][[expUnit]]
                aMat <- matrix(1,nrow=nrow(plotData),ncol=ncol(plotData))
                colnames(aMat) <- colnames(plotData)
                for (gv in gv1) {
                    annLook <- which(plotData[,nc1] >= gv)
                    if (length(annLook)) {
                        aMat[annLook,nc1] <- (aMat[annLook,nc1] + 1)
                    }
                }
                look1 <- which(aMat[,nc1] == as.numeric(ex1))
                if (nrow(plotData)) {
                    colnames(plotData) <- c("x")
                    if (ex1 == "1") pex1 <- "- "
                    else if (ex1 == ux1) pex1 <- "+ "
                    else pex1 <- paste0("_",ex1,"_",ux1," ")
                    p <- ggplot(plotData,aes(x=x)) + geom_histogram(bins=30) + theme_bw(base_size=16)+
                        xlab(nc1)+
                        geom_vline(xintercept=gv1,linetype="dashed",color="red") +
                        theme(legend.position="none")+
                        coord_cartesian(xlim=c(c1qs[1],c1qs[2]))+
                        ggtitle(paste0(nc1,pex1,": ",
                                       length(look1),"/",nrow(plotData)))
                    plotList <- append(plotList,list(p))
                }
            }
        }
        pj <- cowplot::plot_grid(plotlist=plotList)# + theme_cowplot(font_size=24)
        title <- cowplot::ggdraw() + cowplot::draw_label(paste0("Sample ",sampleName), fontface='bold')
        pjOut <- cowplot::plot_grid(title, pj, ncol=1, rel_heights=c(0.1, 1))
        return(pjOut)
    }
    else {
        print("Sample not found in the analysis map.")
        return(NA)
    }
}
