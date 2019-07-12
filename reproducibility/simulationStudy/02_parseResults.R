library(cowplot)
library(viridis)
library(dplyr)
library(tidyr)
library(ggplot2)
#
#Simulation results are written to a local directory. Read in whatever is available.
#
experiments <- list.files("./simResults")
allResults <- matrix(nrow=0,ncol=6)
colnames(allResults) <- c("rate_0.55","rate_0.60","rate_0.65","rate_0.70","rate_0.75","rate_0.80")
faustStatistics <- matrix(nrow=0,ncol=5)
colnames(faustStatistics) <- c("Median # Clusters","True Clusters Found","Median Correlation",
                               "Spiked Correlation","spikeFound")
minFDRCount <- function(pvalVec) {
    qps <- p.adjust(pvalVec,"BH")
    return(length(which(qps==min(qps))))
}
#
#Specify the possible simulation parameters.
#
possibleClusters <- c(10,25,50)
noiseDim <- c(0,1)
batchEffect <- c(0,1)
transType <- c("Identity","Square","Gamma")

getMedBestResults <- function(resultsMatrix) {
    #
    #this function is highly unportable. depends on exact ordering of results from simulation parsing.
    #
    tmp <- matrix(nrow=nrow(resultsMatrix),ncol=0)
    for (i in seq(3,ncol(resultsMatrix),by=4)) {
        tmp <- cbind(tmp,apply(resultsMatrix[,c(i,(i+1))],1,min))
    }
    return(apply(tmp,1,function(x){median(x,na.rm=TRUE)}))
}


firstExp <- TRUE
for (expSetting in experiments){
    if (file.exists(paste0("./simResults/",expSetting,"/intermediateSimOutput.rds"))){
    simRes <- readRDS(paste0("./simResults/",expSetting,"/intermediateSimOutput.rds"))
    expStr <- gsub("c1","10 Clusters\n",expSetting)
    expStr <- gsub("c2","25 Clusters\n",expStr)
    expStr <- gsub("c3","50 Clusters\n",expStr)
    expStr <- gsub("b1","Batch effect: No\n",expStr)
    expStr <- gsub("b2","Batch effect: Yes\n",expStr)
    expStr <- gsub("n1","Nuisance Variables: No\n",expStr)
    expStr <- gsub("n2","Nuisance Variables: Yes\n",expStr)
    expStr <- gsub("t1","Identity Map",expStr)
    expStr <- gsub("t2","Square Map",expStr)
    expStr <- gsub("t3","Gamma Map",expStr)
    clusterNum <- possibleClusters[as.numeric(substr(expSetting,2,2))]
    usedNoise <- noiseDim[as.numeric(substr(expSetting,4,4))]
    usedBatch <- noiseDim[as.numeric(substr(expSetting,6,6))]
    transformationType <- transType[as.numeric(substr(expSetting,8,8))]
    opList <- flowsomList <- resultsList <- list()
    resultsMatrix <- matrix(nrow=0,ncol=7)
    colnames(resultsMatrix) <- c("pctFAUSTUnclass","varSelection","faustTotalCluster",
                                 "faustFoundSpike","spikedCorrelation",
                                 "faustTruePops","trueCorrelation")

    clusterScoresMatrix <- matrix(nrow=14,ncol=0)
    rownames(clusterScoresMatrix) <- c("Purity","Entropy","NormMutInfo",
                                       "VarInfo","NormVarInfo","Specificity",
                                       "Sensitivity","Precision","Recall",
                                       "F-measure","Rand Index","Adj Rand Index",
                                       "Jaccard Index","Fowlkes-Mallows Index")
    for (simIter in seq(length(simRes))) {
        currentSim <- simRes[[simIter]]
        names(currentSim[["flowsomModeling"]])
        resultsMatrix <- rbind(resultsMatrix,
                               c(currentSim[["medianFaustPctUnclassified"]],
                                 currentSim[["faustVarSelectionCount"]],
                                 currentSim[["allFaustPops"]],
                                 currentSim[["spikedInFound"]],
                                 currentSim[["spikedCorrelation"]],
                                 currentSim[["truePopsFound"]],
                                 currentSim[["medFaustCorrelation"]]))
        clusterScoresMatrix <- cbind(clusterScoresMatrix,currentSim[["clusterScores"]])
        faustModel <- currentSim[["faustModeling"]]
        results <- Reduce(rbind,lapply(faustModel,function(x){cbind(x[["responseRate"]],
                                                                    (x[["faustMaxName"]]/50),
                                                                    x[["medianSpikedFDR"]])}))
        resultsList <- append(resultsList,list(results))
        flowsomModel <- currentSim[["flowsomModeling"]]
        resultsFlowsom <- Reduce(rbind,lapply(flowsomModel,function(x){cbind(x[["responseRate"]],
                                                                             x[["medFDRflowSOM"]],
                                                                             x[["medOracleMaxNameFDR"]],
                                                                             x[["medOracleMaxPropFDR"]],
                                                                             x[["medOracleMaxBestFDR"]])}))
        
        flowsomList <- append(flowsomList,list(resultsFlowsom))
        overpartitionedModel <- currentSim[["overpartitionedModeling"]]
        resultsOverpartitioned <- Reduce(rbind,lapply(overpartitionedModel,function(x){cbind(x[["responseRate"]],
                                                                                             x[["medFDRoverpartitioned"]],
                                                                                             x[["medOPMaxNameFDR"]],
                                                                                             x[["medOpMaxPropFDR"]],
                                                                                             x[["medOpMaxBestFDR"]])}))
        opList <- append(opList,list(resultsOverpartitioned))
    }    
    #
    #Record overall results
    #
    fMeasures <- c()
    for (cNum in seq(6)) {
        fMeasures <- append(fMeasures,median(clusterScoresMatrix[10,seq(cNum,ncol(clusterScoresMatrix),by=6),drop=TRUE],na.rm=TRUE))
    }
    names(fMeasures) <- paste0("medFmeasure_",colnames(clusterScoresMatrix)[1:6])
    medianStatistics <- apply(resultsMatrix,2,function(x){median(x,na.rm=TRUE)})
    iterStatistics <- append(medianStatistics,fMeasures)
    expStats <- as.data.frame(t(data.frame(iterStatistics)))
    expStats$clusterNum <- clusterNum
    expStats$usedNoise <- usedNoise
    expStats$usedBatch <- usedBatch
    expStats$transformation <- transformationType
    expStats$setting <- expStr
    print(paste0(expStr,": ",length(simRes)))
    #
    #record faust median results.
    #
    faustResMat <- Reduce(cbind,resultsList)
    medianModeling <- data.frame(responseRate=faustResMat[,1],
                                 medPropSpikeBest=apply(faustResMat[,seq(2,ncol(faustResMat),by=3),drop=FALSE],1,
                                                        function(x){median(x,na.rm=TRUE)}),
                                 medFDRSpike=apply(faustResMat[,seq(3,ncol(faustResMat),by=3),drop=FALSE],1,
                                                   function(x){median(x,na.rm=TRUE)}))
    medianModeling$clusterNum <- clusterNum
    medianModeling$usedNoise <- usedNoise
    medianModeling$usedBatch <- usedBatch
    medianModeling$transformation <- transformationType
    medianModeling$method <- "FAUST"
    medianModeling$setting <- expStr
    medianModeling$numIter <- length(simRes)
    #
    #record oracle flowsom results
    #
    flowsomResMat <- Reduce(cbind,flowsomList)
    medianFlowsom <- data.frame(responseRate=flowsomResMat[,1],
                                medFDRflowSOM=apply(flowsomResMat[,seq(2,ncol(flowsomResMat),by=5),drop=FALSE],1,
                                                    function(x){median(x,na.rm=TRUE)}),
                                maxCountFDRflowSOM=apply(flowsomResMat[,seq(3,ncol(flowsomResMat),by=5),drop=FALSE],1,
                                                    function(x){median(x,na.rm=TRUE)}),
                                maxPropFDRflowSOM=apply(flowsomResMat[,seq(4,ncol(flowsomResMat),by=5),drop=FALSE],1,
                                                        function(x){median(x,na.rm=TRUE)}),
                                bestMedFDRflowSOM=apply(flowsomResMat[,seq(5,ncol(flowsomResMat),by=5),drop=FALSE],1,
                                                        function(x){median(x,na.rm=TRUE)}))
    medianFlowsom$clusterNum <- clusterNum
    medianFlowsom$usedNoise <- usedNoise
    medianFlowsom$usedBatch <- usedBatch
    medianFlowsom$transformation <- transformationType
    medianFlowsom$method <- "flowSOM"
    medianFlowsom$setting <- expStr
    medianFlowsom$numIter <- length(simRes)
    #
    #record overpartitioned flowsom results
    #
    opResMat <- Reduce(cbind,opList)
    medianOp <- data.frame(responseRate=opResMat[,1],
                           medFDRop=apply(opResMat[,seq(2,ncol(opResMat),by=5),drop=FALSE],1,
                                          function(x){median(x,na.rm=TRUE)}),
                           maxCountFDRop=apply(opResMat[,seq(3,ncol(opResMat),by=5),drop=FALSE],1,
                                              function(x){median(x,na.rm=TRUE)}),
                           maxPropFDRop=apply(opResMat[,seq(4,ncol(opResMat),by=5),drop=FALSE],1,
                                              function(x){median(x,na.rm=TRUE)}),
                           bestMedFDRop=apply(opResMat[,seq(5,ncol(opResMat),by=5),drop=FALSE],1,
                                              function(x){median(x,na.rm=TRUE)}))
    medianOp$clusterNum <- clusterNum
    medianOp$usedNoise <- usedNoise
    medianOp$usedBatch <- usedBatch
    medianOp$transformation <- transformationType
    medianOp$method <- "op"
    medianOp$setting <- expStr
    medianOp$numIter <- length(simRes)

    if (firstExp) {
        overallStats <- expStats
        overallFAUST <- medianModeling
        overallFlowsom <- medianFlowsom
        overallOP <- medianOp
        firstExp <- FALSE
    }
    else {
        overallStats <- rbind(overallStats,expStats)
        overallFAUST <- rbind(overallFAUST,medianModeling)
        overallFlowsom <- rbind(overallFlowsom,medianFlowsom)
        overallOP <- rbind(overallOP,medianOp)
    }
    }
}

plotFDR <- function(tTypeStr) {
    subA <- overallFAUST[which(overallFAUST[,"transformation"]==tTypeStr),
                        c("responseRate","setting","medPropSpikeBest")]
    colnames(subA) <- c("responseRate","setting","prop")
    subA$line <- "Median proportion\nof simulation trials\nthat the spiked-in\npopulation amongst\nbest FAUST cluster\nin terms of FDR\nadjusted p-value"
    subB <- overallFAUST[which(overallFAUST[,"transformation"]==tTypeStr),
                         c("responseRate","setting","medFDRSpike")]
    colnames(subB) <- c("responseRate","setting","prop")
    subB$line <- "FAUST"
    subC <- overallFlowsom[which(overallFlowsom[,"transformation"]==tTypeStr),
                           c("responseRate","setting","bestMedFDRflowSOM")]
    colnames(subC) <- c("responseRate","setting","prop")
    subC$line <- "flowSOM\noracle"

    subD <- overallOP[which(overallOP[,"transformation"]==tTypeStr)
                     ,c("responseRate","setting","bestMedFDRop")]
    colnames(subD) <- c("responseRate","setting","prop")
    subD$line <- "flowSOM\noverpartitioned"
    plotDF <- rbind(subB,subC,subD)
    pOut <- ggplot(plotDF,aes(x=responseRate,y=prop,color=line,linetype=line))+
        geom_line(size=1)+
        facet_wrap(~setting)+
        geom_hline(yintercept=0.05,color="#3B528BFF",linetype="dotted")+
        geom_hline(yintercept=0.10,color="#21908CFF",linetype="dotted")+
        geom_hline(yintercept=0.20,color="#5DC863FF",linetype="dotted")+
        scale_color_viridis(discrete=TRUE,begin=0.25,end=0.75,option="B")+
        theme_bw(base_size=16) +
        theme(legend.title=element_blank())+
        xlab("Simulated probability of response when spiked-in population is elevated")+
        ylab("Median observed q-value of spiked-in population")+
        ggtitle("")
    return(pOut)
}
p1 <- plotFDR("Identity")
p2 <- plotFDR("Square")
p3 <- plotFDR("Gamma")

cowplot::ggsave("./simPlots/sim_identityFDRPlot.png",
                p1,
                width=12,
                height=12,
                units="in",
                dpi=300)

cowplot::ggsave("./simPlots/sim_squareFDRPlot.png",
                p2,
                width=12,
                height=12,
                units="in",
                dpi=300)

cowplot::ggsave("./simPlots/sim_gammaFDRPlot.png",
                p3,
                width=12,
                height=12,
                units="in",
                dpi=300)

#
#fmeasure summary plot
#
fmeasureDF <- overallStats[,c(which(grepl("medFmeasure",colnames(overallStats))),17,18)]
plotF <- gather(fmeasureDF,key=MethodLong,value=fmeasure,-c("transformation","setting"))
plotF$allObs <- "FAUST annotated subset"
plotF[which(grepl("_All",plotF$MethodLong)),"allObs"] <- "All Observation"
plotF$Method <- "flowSOM\nOverparitioned"
plotF[which(grepl("flowSOM",plotF$MethodLong)),"Method"] <- "flowSOM\nOracle"
plotF[which(grepl("FAUST",plotF$MethodLong)),"Method"] <- "FAUST"
plotF$setting <- gsub("\\n",", ",plotF$setting)
plotF$setting <- paste0(plotF$setting,"\n")
p1 <- ggplot(plotF[which(plotF$allObs=="All Observation"),],
             aes(x=as.factor(setting),y=fmeasure,color=Method,group=Method))+
    geom_line(size=1,aes(linetype=Method))+
    theme_bw(base_size=16)+
    ylab("")+
    theme(axis.text.x = element_blank())+
    scale_color_viridis(discrete=TRUE,option="B",begin=0.3,end=0.7)+
    geom_hline(yintercept=0.9,linetype="dotted",color="red")+
    ggtitle("All Observations")+
    xlab("")+
    coord_cartesian(ylim=c(0,1))
    
p2 <- ggplot(plotF[-which(plotF$allObs=="All Observation"),],
             aes(x=as.factor(setting),y=fmeasure,color=Method,group=Method))+
    geom_line(size=1,aes(linetype=Method))+
    theme_bw(base_size=16)+
    xlab("Simulation setting")+
    ylab("")+
    theme(axis.text.x = element_blank())+
    scale_color_viridis(discrete=TRUE,option="B",begin=0.3,end=0.7)+
    geom_hline(yintercept=0.9,linetype="dotted",color="red")+
    ggtitle("FAUST Annotated Subset")+
    xlab("")+
    coord_cartesian(ylim=c(0,1))

plotF$subByTrans <- as.factor(paste0(plotF$allObs,"_",plotF$transformation))
p3 <- ggplot(plotF,
             aes(x=as.factor(setting),y=fmeasure,color=Method,group=Method))+
    geom_line(size=1,aes(linetype=Method))+
    theme_bw(base_size=16)+
    facet_wrap(~subByTrans)+
    xlab("Simulation setting")+
    ylab("Median F-measure across simulation iterations")+
    theme(axis.text.x = element_blank())+
    scale_color_viridis(discrete=TRUE,option="B",begin=0.3,end=0.7)+
    geom_hline(yintercept=0.9,linetype="dotted",color="red")+
    ggtitle("FAUST Annotated Subset")+
    coord_cartesian(ylim=c(0,1))
pj <- plot_grid((p1+theme(legend.position="none")),(p2+theme(legend.position="none")),nrow=2,rel_heights=c(0.35,1))
ylbl <- ggdraw() + draw_label("                                                                                                               Median F-measure across simulation iterations",
                              fontface='bold',angle=90,size=16)
plbl <- plot_grid(ylbl,pj,nrow=1,rel_widths=c(0.025,1))
legend <- get_legend(p1+theme(legend.position="bottom"))
pOut <- plot_grid(plbl, legend, ncol=1,rel_heights = c(1,0.1))
cowplot::ggsave("./simPlots/sim_fmeasurePlot.png",
                p3,
                width=12,
                height=12,
                units="in",
                dpi=300)

#
#Correlation statistics plot
#
sub2 <- cbind(overallStats[,c("spikedCorrelation","transformation","setting")],"True counts of \nelevated population and\nFAUST  counts of\nelevated population\nmatched by\nannotation\n\n\n")
colnames(sub2) <- c("Score","transformation","Setting","Measure")
sub3 <- cbind(overallStats[,c("trueCorrelation","transformation","setting")],"Median true counts of \nall simulated populations\nvs FAUST counts\nwith matching annotations\n\n\n")
colnames(sub3) <- c("Score","transformation","Setting","Measure")
plotDF <- rbind(sub2,sub3)
plotDF$Setting <- gsub("\\n",", ",plotDF$Setting)
pCor <- ggplot(plotDF,aes(x=as.factor(Setting),y=Score,color=Measure,group=Measure))+
    geom_line(size=1,aes(linetype=Measure))+
    facet_wrap(~transformation)+
    theme_bw(base_size=16)+
    xlab("Simulation setting")+
    ylab("")+
    geom_hline(yintercept=0.9,linetype="dotted",color="red")+
    ylab("Median correlation across simulation iterations")+
    theme(axis.text.x = element_blank())+
    scale_color_viridis(discrete=TRUE,option="B",begin=0.3,end=0.7)+
    coord_cartesian(ylim=c(0,1))
legend <- get_legend(pCor+theme(legend.position="bottom"))
pOut <- plot_grid((pCor+theme(legend.position="none")), legend, ncol=1,rel_heights = c(1,0.1))
cowplot::ggsave("./simPlots/sim_correlationPlot.png",
                pCor,
                width=12,
                height=12,
                units="in",
                dpi=300)

#
#Number of Clusters Plot
#
sub1 <- cbind(overallStats[,c("clusterNum","transformation","setting")],"True number\nof clusters\n\n")
colnames(sub1) <- c("Score","transformation","Setting","Measure")
sub2 <- cbind(overallStats[,c("faustTotalCluster","transformation","setting")],"Median number of\nall FAUST clusters\n\n")
colnames(sub2) <- c("Score","transformation","Setting","Measure")
sub3 <- cbind(overallStats[,c("faustTruePops","transformation","setting")],"Median number of\ntrue clusters FAUST\nannotated\n\n")
colnames(sub3) <- c("Score","transformation","Setting","Measure")
plotDF <- rbind(sub1,sub2,sub3)#,sub4)
plotDF$Setting <- gsub("\\n",", ",plotDF$Setting)
pNoc <- ggplot(plotDF,aes(x=as.factor(Setting),y=Score,color=Measure,group=Measure))+
    geom_line(size=1,aes(linetype=Measure))+
    theme_bw(base_size=16)+
    facet_wrap(~transformation)+
    xlab("Simulation setting")+
    ylab("Median measure across simulation iterations")+
    theme(axis.text.x = element_blank())+
    scale_color_viridis(discrete=TRUE,option="B",begin=0.3,end=0.7)
legend <- get_legend(pNoc+theme(legend.position="bottom"))
pOut <- plot_grid((pNoc+theme(legend.position="none")), legend, ncol=1,rel_heights = c(1,0.1))
cowplot::ggsave("./simPlots/sim_nocPlot.png",
                pNoc,
                width=12,
                height=12,
                units="in",
                dpi=300)

