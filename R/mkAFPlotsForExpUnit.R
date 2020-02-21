#' Fast Annotation Using Shape-constrained Trees
#'
#' mkAFPlotsForExpUnit produces joyplot of marginal and conditional densities along with their gate locations
#' for the specified experimental unit.
#'
#' @param projectPath An absolute path to a directory on your system. Output from the
#' FAUST pipeline is written to the directory "projectPath/faustData".
#'
#' @param expUnit The experimental unit you want to plot.
#'
#' @param threadNum Number of threads to use.
#'
#' @param debugFlag Boolean value. Set to TRUE to print status of plotting code.
#'
#' @param maxRidgesAtDepth Numeric value. Maximum number of ridgelines to plot at any depth.
#' In very large annotation forest, this can be used to prevent overplotting in derived
#' ridgeline plot.
#'
#' @export
#' @md
#' @importFrom cowplot save_plot
#' @importFrom stats density
#' @importFrom ggplot2 ggplot geom_segment xlab ylab theme theme_bw ggtitle element_blank coord_cartesian aes_string
#' @importFrom ggridges geom_density_ridges
#' @importFrom viridis scale_fill_viridis
#' @importFrom dplyr group_by summarize left_join
#' @importFrom stats density
#' @import tidyr
mkAFPlotsForExpUnit <- function(projectPath,
                                expUnit,
                                threadNum = 1,
                                debugFlag = FALSE,
                                maxRidgesAtDepth=Inf)
{
    if (!dir.exists(file.path(projectPath,"faustData","plotData"))) {
        print(paste0("faustData/plotData directory not detected in ",projectPath))
        print("This indicates the FAUST pipeline has not been run.")
        print("If it has, update the projectPath to the directory containing faustData")
        print("If it has not, run the FAUST pipeline before calling this function.")
        return(NA)
    }
    if (!dir.exists(file.path(normalizePath(projectPath),
                              "faustData",
                              "plotData",
                              "afPlots"))) {
        dir.create(file.path(normalizePath(projectPath),
                             "faustData",
                             "plotData",
                             "afPlots"))
    }
    if (!dir.exists(file.path(normalizePath(projectPath),
                              "faustData",
                              "plotData",
                              "afPlots",
                              expUnit))) {
        dir.create(file.path(normalizePath(projectPath),
                             "faustData",
                             "plotData",
                             "afPlots",
                             expUnit))
    }
    expUnitExprs <- readRDS(file.path(normalizePath(projectPath),
                                    "faustData",
                                    "expUnitData",
                                    expUnit,
                                    "expUnitExprs.rds"))
    expUnitRes <- readRDS(file.path(normalizePath(projectPath),
                                  "faustData",
                                  "expUnitData",
                                  expUnit,
                                  "expUnitRes.rds"))
    channelBounds <- readRDS(file.path(normalizePath(projectPath),
                                       "faustData",
                                       "metaData",
                                       "channelBounds.rds"))
    resFlag <- as.logical(max(apply(expUnitRes,2,max)))
    annF <- growAnnotationForest(dataSet = expUnitExprs,
                                 numberIterations = 1,
                                 pValueThreshold = 0.25,
                                 minimumClusterSize = 25,
                                 randomCandidateSearch = FALSE,
                                 maximumSearchDepth = 2,
                                 numberOfThreads = threadNum,
                                 maximumGatingNum = 1e10,
                                 anyValueRestricted = resFlag,
                                 resValMatrix = expUnitRes,
                                 cutPointUpperBound = 2,
                                 getDebugInfo = debugFlag,
                                 randomSeed = 20180917,
                                 recordCounts = FALSE,
                                 recordIndices = TRUE)
    for (channel in names(annF[["indexData"]])) {
        if (debugFlag) print(paste0("Drawing plot for ",channel))
        outChannel <- gsub("[[:punct:]]","",channel)
        outChannel <- gsub("[[:space:]]","",outChannel)
        outChannel <- gsub("[[:cntrl:]]","",outChannel)
        lowBound <- channelBounds["Low",channel]
        highBound <- channelBounds["High",channel]
        allIndices <- annF[["indexData"]][[channel]]
        allDepths <- annF[["indexDepthData"]][[channel]]
        numObs <- nrow(expUnitExprs)
        startI <- 1
        endI <- numObs
        endLength <- length(allIndices)
        indexLengths <- c()
        expressionDepths <- expressionData <- c()
        indexCtr <- 1
        gateList <- list()
        depthCtr <- c()
        while (startI < endLength) {
            currentLookup <- allIndices[startI:endI]
            allExprs <- expUnitExprs[currentLookup,channel]
            lowLookup <- which(allExprs <= max(lowBound,-Inf))
            highLookup <- which(allExprs >= min(highBound,Inf))
            plotLookup <- setdiff(seq(length(allExprs)),c(sort(unique(lowLookup,highLookup))))
            plotSummary <- table(depthCtr)
            if (length(plotLookup)) {
                curExprs <- allExprs[plotLookup]
                curDepth <- rep(allDepths[indexCtr],length(curExprs))
                if ((length(which(names(plotSummary) == curDepth[1]))==0)||
                    (plotSummary[which(names(plotSummary) == curDepth[1])]<
                     maxRidgesAtDepth)) {
                    depthCtr <- append(depthCtr,curDepth[1])
                    expressionData <- append(expressionData,curExprs)
                    expressionDepths <- append(expressionDepths,curDepth)
                    gateLocs <- tsGates(sort(addNoiseToDataVector(curExprs,4,123)),0)
                    gateLocs <- gateLocs[-c(1,length(gateLocs))]
                    gateList <- append(gateList,list(gateLocs))
                    names(gateList)[length(gateList)] <- indexCtr
                    indexLengths <- append(indexLengths,rep(indexCtr,length(plotLookup)))
                    indexCtr <- indexCtr + 1
                }
            }
            startI <- endI + 1
            endI <- startI + numObs - 1
        }
        #add the margin
        allExprs <- expUnitExprs[,channel]
        lowLookup <- which(allExprs <= max(lowBound,-Inf))
        highLookup <- which(allExprs >= min(highBound,Inf))
        plotLookup <- setdiff(seq(length(allExprs)),c(sort(unique(lowLookup,highLookup))))
        curExprs <- allExprs[plotLookup]
        expressionData <- append(expressionData,curExprs)
        curDepth <- rep(1,length(curExprs))
        expressionDepths <- append(expressionDepths,curDepth)
        gateLocs <- -Inf
        gateList <- append(gateList,list(gateLocs))
        names(gateList)[length(gateList)] <- indexCtr
        indexLengths <- append(indexLengths,rep(indexCtr,length(plotLookup)))
        indexCtr <- indexCtr + 1
        maxLen <- max(unlist(lapply(gateList,length)))
        for (gateNum in seq(length(gateList))) {
            cLen <- length(gateList[[gateNum]])
            if (cLen < maxLen)
                gateList[[gateNum]] <- c(gateList[[gateNum]],rep(gateList[[gateNum]][1], maxLen - cLen))
        }
        gateDF <- as.data.frame(Reduce(rbind,gateList))
        colnames(gateDF) <- paste0("x",seq(ncol(gateDF)))
        keepNames <- colnames(gateDF)
        #two-line to remove note in R CMD check
        yData <- names(gateList)
        gateDF$yData <- yData
        if (maxRidgesAtDepth < Inf) {
            #only save these files if we are limiting the number of ridges
            #assumes user wishes to generate plots with a requested sub-collection
            saveRDS(indexCtr,
                    file.path(normalizePath(projectPath),
                              "faustData",
                              "plotData",
                              "afPlots",
                              expUnit,
                              paste0(outChannel,"_indexCtr.rds")))
            saveRDS(indexLengths,
                    file.path(normalizePath(projectPath),
                              "faustData",
                              "plotData",
                              "afPlots",
                              expUnit,
                              paste0(outChannel,"_indexLengths.rds")))
            saveRDS(expressionData,
                    file.path(normalizePath(projectPath),
                              "faustData",
                              "plotData",
                              "afPlots",
                              expUnit,
                              paste0(outChannel,"_expressionData.rds")))
            saveRDS(expressionDepths,
                    file.path(normalizePath(projectPath),
                              "faustData",
                              "plotData",
                              "afPlots",
                              expUnit,
                              paste0(outChannel,"_expressionDepths.rds")))
            saveRDS(gateDF,
                    file.path(normalizePath(projectPath),
                              "faustData",
                              "plotData",
                              "afPlots",
                              expUnit,
                              paste0(outChannel,"_gateDF.rds")))
        }
        indexScaling <- c()
        for (indexNum in seq(indexCtr - 1)) {
            indexScaling <- append(indexScaling,length(which(indexLengths == indexNum)))
            names(indexScaling)[length(indexScaling)] <- indexNum
        }
        indexScaling <- indexScaling/min(indexScaling)
        for (indexNum in seq(indexCtr - 1)) {
            scaleVal <- as.numeric(indexScaling[as.character(indexNum)])
            densityEstimate <- density(expressionData[which(indexLengths == indexNum)],bw = "SJ")
            xData <- densityEstimate$x
            hData <- (densityEstimate$y*scaleVal)
            newDF <- data.frame(xData = xData,
                                yData = as.character(indexNum),
                                hData = hData,
                                stringsAsFactors = FALSE)
            if (indexNum == 1) {
                plotDF <- newDF
            }
            else {
                plotDF <- rbind(plotDF,newDF)
            }
        }
        mhDF <- as.data.frame(group_by(plotDF,yData) %>% summarize(meanHeight = mean(hData)) )
        mhDF2 <- mhDF[order(mhDF[,"meanHeight"]),]
        #assing y1 to NULL for R CMD check note
        y1 <- NULL
        mhDF2$y1 <- seq(nrow(mhDF2))
        mhDF2$alphaScale <- seq(nrow(mhDF2))/nrow(mhDF2)
        plotDF2 <- left_join(plotDF,mhDF2,by="yData")
        gateDF2 <- left_join(gateDF,mhDF2,by="yData")
        gateDF2 <- gateDF2[,c(keepNames,"y1")]
        gateDF3 <- gather(data=gateDF2,key="xName",value="x1",-y1)
        gateDF3$x2 <- gateDF3$x1
        gateDF3$y2 <- as.numeric(gateDF3$y1 + 1)
        saveRDS(plotDF2,
                file.path(normalizePath(projectPath),
                          "faustData",
                          "plotData",
                          "afPlots",
                          expUnit,
                          paste0(outChannel,"_plotDF2.rds")))
        saveRDS(gateDF3,
                file.path(normalizePath(projectPath),
                          "faustData",
                          "plotData",
                          "afPlots",
                          expUnit,
                          paste0(outChannel,"_gateDF3.rds")))
        p <- ggplot(plotDF2,aes_string(x = "xData", y = "y1", height = "hData",
                                       group = "y1", fill = "y1"))+
            geom_density_ridges(data = plotDF2,
                                scale = 10,
                                alpha = plotDF2$alphaScale,
                                stat = "identity",
                                rel_min_height = 0.001) +
            geom_segment(data = gateDF3,
                         aes_string(x = "x1", xend = "x2", y = "y1", yend = "y2"),
                         color = "red", inherit.aes = FALSE) +
            xlab("") +
            scale_fill_viridis() +
            ylab("") +
            theme_bw() +
            theme(axis.ticks.y = element_blank(),
                  axis.text.y = element_blank(),
                  legend.position = "none")
        pOut <- p + ggtitle(paste0(channel,": Annotation Forest with Gates")) + xlab("Expression Value")
        save_plot(file.path(normalizePath(projectPath),
                            "faustData",
                            "plotData",
                            "afPlots",
                            expUnit,
                            paste0(outChannel,".png")),
                  pOut,
                  base_width = 15,
                  base_height = 15)
    }
    return()
}
