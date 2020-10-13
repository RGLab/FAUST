#' Full Annotation Using Shape-constrained Trees
#'
#' backwardPhenotypeSelection
#' 
#' @param projectPath A path to a directory containing a completed FAUST analysis. 
#' The projectPath directory is expected to have a faustData sub-directory, which
#' contains a completed faust analysis. 
#'
#' @param startingPhenotype A character vector containing the phenotype of interest.
#' This must match a column in the file 'faustCountMatrix.rds' that is produced
#' by a complete faust run.
#'
#' @param inputModelStr A string containing a model formula. This formula is used
#' to guide the backward selection. The LHS of the formula must read
#' 'cbind(childCount,(parentCount-childCount))'. The RHS of the formula
#' (after the '~' character) must match column names of the metaDataDF.
#' The variable of interest must occupy the first position of the 
#'
#' @param metaDataDF A data frame with all relevant meta data (except the counts)
#' needed to fit the model formula described in the 'inputModelStr' parameter.
#' One variable must be called 'sampleName', and the entries must match folders
#' contained in projectPath/faustData/sampleData.
#' 
#' @return A list with two slots.
#'
#' Slot 'selections': the markers chosen at each selection stage. 
#'
#' Slot 'pvalues': the associated pvalues.
#' 
#' @export
#' @md
#' @importFrom lme4 glmer
#' @importFrom dplyr inner_join
backwardPhenotypeSelection <- function(
                                       projectPath,
                                       startingPhenotype,
                                       inputModelStr,
                                       metaDataDF
                                       )
{
    #
    #set up the environment
    #
    encodingDF <- readRDS(file.path(projectPath,"faustData","metaData","exhaustiveColNameMap.rds"))
    internalPhenotype <- encodingDF[which(encodingDF[,"newColNames"]==startingPhenotype),"faustColNames"]
    phenoData <- strsplit(internalPhenotype,"~")[[1]]
    possibleMarkers <- phenoData[seq(1,length(phenoData),by=3)]
    #
    #targetExrps will contain all possible markers to target
    #and their associated expression level for looking up in the annotation
    #matrix
    #
    targetExprs <- phenoData[seq(2,length(phenoData),by=3)]
    names(targetExprs) <- possibleMarkers
    #
    #get the markers in the column order of the annotation matrix.
    #
    gatingFiles <- list.files(file.path(projectPath,"faustData","gateData"))
    selCFN <- gatingFiles[grepl("_selectedChannels.rds",gatingFiles)]
    selectedChannels <- readRDS(file.path(projectPath,"faustData","gateData",selCFN))
    #
    #first, get the p-value for the complete phenotype
    #
    admissibleMarkers <- selectedChannels
    phenoOutList <- bwMarkerOutList <- pvalOutList <- list()
    modelPval <- .getPvalueForMarkers(
        projectPath=projectPath,
        candidateMarkers=admissibleMarkers,
        targetExprs=targetExprs,
        selectedChannels=selectedChannels,
        metaDataDF=metaDataDF,
        modelStr=inputModelStr
    )
    #
    #record the data for the complete phenotype
    #
    bwMarkerOutList <- append(bwMarkerOutList,list(admissibleMarkers))
    names(bwMarkerOutList)[length(bwMarkerOutList)] <- paste0(length(admissibleMarkers)," markers")
    pvalOutList <- append(pvalOutList,list(modelPval))
    names(pvalOutList)[length(pvalOutList)] <- paste0(length(admissibleMarkers)," markers")
    phenoOutList <- append(phenoOutList,list(startingPhenotype))
    names(phenoOutList)[length(phenoOutList)] <- paste0(length(admissibleMarkers)," markers")
    #
    #iterate over the markers, eliminating markers in a step-wise fashion.
    #
    while (length(admissibleMarkers) > 1) {
        print(paste0("Processing ",length(admissibleMarkers)," markers"))
        pvalVec <- c()
        markerList <- list()
        for (cmn in admissibleMarkers) {
            candMarkers <- setdiff(admissibleMarkers,cmn)
            modelPval <- .getPvalueForMarkers(
                projectPath=projectPath,
                candidateMarkers=candMarkers,
                targetExprs=targetExprs,
                selectedChannels=selectedChannels,
                metaDataDF=metaDataDF,
                modelStr=inputModelStr
            )
            pvalVec <- append(pvalVec,modelPval)
            markerList <- append(markerList,list(candMarkers))
        }
        #
        #select the choice of markers that minimizes the pvalues
        #in case of ties, select the first choice.
        #
        cbPvalue <- pvalVec[[which(pvalVec==min(pvalVec))[1]]]
        cbSelection <- markerList[[which(pvalVec==min(pvalVec))[1]]]
        #
        #get the expression associated with the selected markers
        #
        markerVec <- c()
        for (smn in cbSelection) {
            markerIndex <- which(grepl(smn,phenoData))
            markerStr <- paste0(paste0(phenoData[seq(markerIndex,(markerIndex+2))],collapse="~"),"~")
            markerStr <- gsub("~1~2~","- ",markerStr)
            markerStr <- gsub("~2~2~","+ ",markerStr)
            markerStr <- gsub("~1~3~","Low ",markerStr)
            markerStr <- gsub("~2~3~","Dim ",markerStr)
            markerStr <- gsub("~3~3~","Bright ",markerStr)
            markerVec <- append(markerVec,markerStr)
        }
        markerVecOut <- paste0(markerVec,collapse="")
        #
        #record this round of selection
        #
        bwMarkerOutList <- append(bwMarkerOutList,list(cbSelection))
        names(bwMarkerOutList)[length(bwMarkerOutList)] <- paste0(length(cbSelection)," markers")
        pvalOutList <- append(pvalOutList,list(cbPvalue))
        names(pvalOutList)[length(pvalOutList)] <- paste0(length(cbSelection)," markers")
        phenoOutList <- append(phenoOutList,list(markerVecOut))
        names(phenoOutList)[length(phenoOutList)] <- paste0(length(cbSelection)," markers")
        #
        #update the admissible candidates
        #
        admissibleMarkers <- cbSelection
    }
    #
    #store the results and return.
    #
    outList <- list(
        `selections`=bwMarkerOutList,
        `phenotypes`=phenoOutList,
        `pvalues`=pvalOutList
    )
    return(outList)
}

.modelContrast <- function(modelStr,dataSet) {
    out <- tryCatch(
    {
        m <- glmer(
            formula=as.formula(modelStr),
            data=dataSet,
            family="binomial",
            control=glmerControl(
                optimizer="bobyqa", optCtrl = list(maxfun = 1e9),
                check.conv.singular = .makeCC(action = "warning",  tol = 1e-4))
        )
        rv <- c(coefficients(summary(m))[2,4])
        return(rv)
    },
    error=function(cond){
        message("Error!")
        message(cond)
        return(NA)
    },
    warning=function(cond){
        message("Warning!")
        message(cond)
        return(NA)
    },
    finally={
    }
    )
    return(out)
}

.getPhenotypeCountDF <- function(projectPath,
                                 markersInPheno,
                                 targetExprs,
                                 selectedChannels)
{
    allSampleNames <- list.files(file.path(projectPath,"faustData","sampleData"))
    pcVec <- ccVec <- snVec <- c()
    for (sampleName in allSampleNames) {
        #
        #load the per-cell annotation matrix
        #
        amIn <- read.table(
            file=file.path(projectPath,"faustData","sampleData",sampleName,"annotationMatrix.csv"),
            header=FALSE,
            sep=",",
            stringsAsFactors=FALSE
        )
        colnames(amIn) <- selectedChannels
        parentCount <- nrow(amIn)
        #
        #get counts for the markers that constitute the phenotype
        #
        childLookup <- seq(nrow(amIn))
        for (mn in markersInPheno) {
            currentLookup <- which(amIn[,mn]==as.numeric(targetExprs[[mn]]))
            childLookup <- intersect(childLookup,currentLookup)
        }
        childCount <- length(childLookup)
        #
        #record data from sample
        #
        pcVec <- append(pcVec,parentCount)
        ccVec <- append(ccVec,childCount)
        snVec <- append(snVec,sampleName)
    }
    phenotypeCountDF <- data.frame(
        childCount=ccVec,
        parentCount=pcVec,
        sampleName=snVec,
        stringsAsFactors=FALSE
    )
    return(phenotypeCountDF)
}

.getPvalueForMarkers <- function(projectPath,
                                 candidateMarkers,
                                 targetExprs,
                                 selectedChannels,
                                 metaDataDF,
                                 modelStr)
{
    #
    #get counts for the phenotype consisting of the candidate markers
    #
    phenoDF <- .getPhenotypeCountDF(projectPath=projectPath,
                                    markersInPheno=candidateMarkers,
                                    targetExprs=targetExprs,
                                    selectedChannels=selectedChannels)
    #
    #merge on the meta data, fit the model, and resut
    #
    modelDF <- inner_join(phenoDF,metaDataDF,by=c("sampleName"))
    modelP <- .modelContrast(
        modelStr=modelStr,
        dataSet=modelDF
    )
    return(modelP)
}
