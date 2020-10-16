#' backwardPhenotypeSelection
#' 
#' @param projectPath A path to a directory containing a completed FAUST analysis. 
#' The projectPath directory is expected to have a faustData sub-directory, which
#' contains a completed faust analysis. 
#'
#' @param startingPhenotype A character vector specifying the phenotype of interest.
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
#' @param markersInParentPhenotype The subset of markers in the startingPhenotype
#' used to compute the sample-by-sample parentCount. NULL by default, leads
#' to using all observations in each sample as the parentCount. Example:
#' if the startingPhenotype contains "CD3+ CD4+ CD8-", setting this parameter to
#' c("CD3","CD4","CD8") will cause the parent count in the GLMM to be total
#' CD4 T cells rather than all cells in the sample.
#' 
#' @return A list with four slots.
#'
#' Slot 'selections': the markers chosen at each selection stage. 
#'
#' Slot 'phenotypes': the phenotype associated with the selected markers.
#'
#' Slot 'pvalues': the associated pvalues.
#'
#' Slot `encodingData`: data used to construct phenotype queries based on the
#' results.
#'
#' @importFrom lme4 glmer
#' @importFrom dplyr inner_join
#' @importFrom bit bitwhich
#' @export
backwardPhenotypeSelection <- function(
                                       projectPath,
                                       startingPhenotype,
                                       inputModelStr,
                                       metaDataDF,
                                       markersInParentPhenotype=NULL
                                       )
{
    require(bit)  #for the bitwhich class, and associated operations
    require(lme4) #for glmer
    require(dplyr) #for inner_join
    #
    #set up the environment
    #
    encodingDF <- readRDS(file.path(projectPath,"faustData","metaData","exhaustiveColNameMap.rds"))
    internalPhenotype <- encodingDF[which(encodingDF[,"newColNames"]==startingPhenotype),"faustColNames"]
    phenoData <- strsplit(internalPhenotype,"~")[[1]]
    possibleMarkers <- phenoData[seq(1,length(phenoData),by=3)]
    #
    #test to make sure markers in the parent phenotype are a subset of the startingPhenotype
    #
    if (!is.null(markersInParentPhenotype)) {
        for (pmn in markersInParentPhenotype) {
            if (length(which(grepl(paste0("\\<",pmn,"\\>"),phenoData)))==0) {
                stop("Phenotypes in markersInParentPhenotype must be a subset of those in the startingPhenotype")
            }
        }
    }
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
    #next, get the phenotype lookups across the dataaset
    #
    phenotypeLookupList <- .getPhenotypeLookupList(
        projectPath=projectPath,
        targetExprs=targetExprs,
        selectedChannels=selectedChannels
    )
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
        modelStr=inputModelStr,
        phenotypeLookupList=phenotypeLookupList,
        markersInParentPhenotype=markersInParentPhenotype
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
    while (length(admissibleMarkers) > (length(markersInParentPhenotype)+1)) { #test relies on length(NULL)==0
        print(paste0("Processing ",length(admissibleMarkers)," markers"))
        pvalVec <- c()
        markerList <- list()
        for (cmn in setdiff(admissibleMarkers,markersInParentPhenotype)) {
            candMarkers <- setdiff(admissibleMarkers,cmn)
            modelPval <- .getPvalueForMarkers(
                projectPath=projectPath,
                candidateMarkers=candMarkers,
                targetExprs=targetExprs,
                selectedChannels=selectedChannels,
                metaDataDF=metaDataDF,
                modelStr=inputModelStr,
                phenotypeLookupList=phenotypeLookupList ,
                markersInParentPhenotype=markersInParentPhenotype
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
            markerIndex <- which(grepl(paste0("\\<",smn,"\\>"),phenoData))
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
        `pvalues`=pvalOutList,
        `encodingData`=phenoData
    )
    return(outList)
}


#' forwardPhenotypeSelection
#' 
#' @param projectPath A path to a directory containing a completed FAUST analysis. 
#' The projectPath directory is expected to have a faustData sub-directory, which
#' contains a completed faust analysis. 
#'
#' @param targetPhenotype A character vector specifying the phenotype of interest.
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
#' @param markersInParentPhenotype The subset of markers in the targetPhenotype 
#' used to compute the sample-by-sample parentCount. NULL by default, leads
#' to using all observations in each sample as the parentCount. Example:
#' if the targetPhenotype contains "CD3+ CD4+ CD8-", setting this parameter to
#' c("CD3","CD4","CD8") will cause the parent count in the GLMM to be total
#' CD4 T cells rather than all cells in the sample.
#' 
#' @return A list with four slots.
#'
#' Slot 'selections': the markers chosen at each selection stage. 
#'
#' Slot 'phenotypes': the phenotype associated with the selected markers.
#'
#' Slot 'pvalues': the associated pvalues.
#'
#' Slot `encodingData`: data used to construct phenotype queries based on the
#' results.
#'
#' @importFrom lme4 glmer
#' @importFrom dplyr inner_join
#' @importFrom bit bitwhich
#' @export
forwardPhenotypeSelection <- function(
                                      projectPath,
                                      targetPhenotype,
                                      inputModelStr,
                                      metaDataDF,
                                      markersInParentPhenotype=NULL
                                      )
{
    require(bit)  #for the bitwhich class, and associated operations
    require(lme4) #for glmer
    require(dplyr) #for inner_join
    #
    #set up the environment
    #
    encodingDF <- readRDS(file.path(projectPath,"faustData","metaData","exhaustiveColNameMap.rds"))
    internalPhenotype <- encodingDF[which(encodingDF[,"newColNames"]==targetPhenotype),"faustColNames"]
    phenoData <- strsplit(internalPhenotype,"~")[[1]]
    possibleMarkers <- phenoData[seq(1,length(phenoData),by=3)]
    #
    #test to make sure markers in the parent phenotype are a subset of the targetPhenotype
    #
    if (!is.null(markersInParentPhenotype)) {
        for (pmn in markersInParentPhenotype) {
            if (length(which(grepl(paste0("\\<",pmn,"\\>"),phenoData)))==0) {
                stop("Phenotypes in markersInParentPhenotype must be a subset of those in the targetPhenotype")
            }
        }
    }
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
    #next, get the phenotype lookups across the dataaset
    #
    phenotypeLookupList <- .getPhenotypeLookupList(
        projectPath=projectPath,
        targetExprs=targetExprs,
        selectedChannels=selectedChannels
    )
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
        modelStr=inputModelStr,
        phenotypeLookupList=phenotypeLookupList,
        markersInParentPhenotype=markersInParentPhenotype
    )
    #
    #record the data for the complete phenotype
    #
    bwMarkerOutList <- append(bwMarkerOutList,list(admissibleMarkers))
    names(bwMarkerOutList)[length(bwMarkerOutList)] <- paste0(length(admissibleMarkers)," markers")
    pvalOutList <- append(pvalOutList,list(modelPval))
    names(pvalOutList)[length(pvalOutList)] <- paste0(length(admissibleMarkers)," markers")
    phenoOutList <- append(phenoOutList,list(targetPhenotype))
    names(phenoOutList)[length(phenoOutList)] <- paste0(length(admissibleMarkers)," markers")
    #
    #iterate over the markers, adding markers in a step-wise fashion.
    #    
    if (is.null(markersInParentPhenotype)) {
        stepwiseSelection <- c()
    }
    else {
        stepwiseSelection <- markersInParentPhenotype
        admissibleMarkers <- setdiff(selectedChannels,markersInParentPhenotype)
    }
    while (length(admissibleMarkers) > 1) {
        print(paste0("Processing ",(length(stepwiseSelection)+1)," markers"))
        pvalVec <- c()
        markerList <- list()
        for (cmn in admissibleMarkers) {
            candMarkers <- append(stepwiseSelection,cmn)
            modelPval <- .getPvalueForMarkers(
                projectPath=projectPath,
                candidateMarkers=candMarkers,
                targetExprs=targetExprs,
                selectedChannels=selectedChannels,
                metaDataDF=metaDataDF,
                modelStr=inputModelStr,
                phenotypeLookupList=phenotypeLookupList,
                markersInParentPhenotype=markersInParentPhenotype
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
            markerIndex <- which(grepl(paste0("\\<",smn,"\\>"),phenoData))
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
        #record selection path and update the admissible candidates
        #
        stepwiseSelection <- cbSelection
        admissibleMarkers <- setdiff(selectedChannels,stepwiseSelection)
    }
    #
    #store the results and return.
    #
    outList <- list(
        `selections`=bwMarkerOutList,
        `phenotypes`=phenoOutList,
        `pvalues`=pvalOutList,
        `encodingData`=phenoData
    )
    return(outList)
}



#' getCountsForTargetMarkers
#' 
#' @param projectPath A path to a directory containing a completed FAUST analysis. 
#' The projectPath directory is expected to have a faustData sub-directory, which
#' contains a completed faust analysis. 
#'
#' @param referencePhenotype A character vector specifying the phenotype of interest.
#' The expression of this phenotype is used by the target markers to derive the count 
#' matrix. This phenotype must match a column in the file 'faustCountMatrix.rds' 
#' that is produced by a complete faust run.
#'
#' @param targetMarkers The markers of interest in the sub-phenotype.
#'
#' @param metaDataDF A data frame with all relevant meta data (except the counts)
#' needed to fit the model formula described in the 'inputModelStr' parameter.
#' One variable must be called 'sampleName', and the entries must match folders
#' contained in projectPath/faustData/sampleData.
#'
#' @param markersInParentPhenotype The subset of markers in the referencePhenotype
#' used to compute the sample-by-sample parentCount. NULL by default, leads
#' to using all observations in each sample as the parentCount. Example:
#' if the referencePhenotype contains "CD3+ CD4+ CD8-", setting this parameter to
#' c("CD3","CD4","CD8") will cause the parent count in the GLMM to be total
#' CD4 T cells rather than all cells in the sample.
#' 
#' @return A data frame with counts for the sub-phenotype consisting of the target
#' markers. All meta data in the metaDataDF is merged on.
#'
#' @importFrom dplyr inner_join
#' @importFrom bit bitwhich
#' @export
getCountsForTargetMarkers <- function(
                                      projectPath,
                                      referencePhenotype,
                                      targetMarkers,
                                      metaDataDF,
                                      markersInParentPhenotype=NULL
                                      )
{
    require(bit)  #for the bitwhich class, and associated operations
    require(lme4) #for glmer
    require(dplyr) #for inner_join
    #
    #set up the environment
    #
    encodingDF <- readRDS(file.path(projectPath,"faustData","metaData","exhaustiveColNameMap.rds"))
    internalPhenotype <- encodingDF[which(encodingDF[,"newColNames"]==referencePhenotype),"faustColNames"]
    phenoData <- strsplit(internalPhenotype,"~")[[1]]
    possibleMarkers <- phenoData[seq(1,length(phenoData),by=3)]
    #
    #test to make sure markers in the parent phenotype are a subset of the referencePhenotype
    #
    if (!is.null(markersInParentPhenotype)) {
        for (pmn in markersInParentPhenotype) {
            if (length(which(grepl(paste0("\\<",pmn,"\\>"),phenoData)))==0) {
                stop("Phenotypes in markersInParentPhenotype must be a subset of those in the referencePhenotype")
            }
        }
    }
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
    #next, get the phenotype lookups across the dataaset
    #
    phenotypeLookupList <- .getPhenotypeLookupList(
        projectPath=projectPath,
        targetExprs=targetExprs,
        selectedChannels=selectedChannels
    )
    #
    #get counts for the sub-phenotype consisting of the candidate markers
    #
    phenoDF <- .getPhenotypeCountDF(
        projectPath=projectPath,
        markersInPheno=targetMarkers,
        targetExprs=targetExprs,
        selectedChannels=selectedChannels,
        phenotypeLookupList=phenotypeLookupList,
        markersInParentPhenotype=markersInParentPhenotype
    )
    #
    #merge on the meta data, and return the data frame
    #
    outputDF <- inner_join(phenoDF,metaDataDF,by=c("sampleName"))
    return(outputDF)
}


.getPhenotypeLookupList <- function(projectPath,
                                   targetExprs,
                                   selectedChannels)
{
    allSampleNames <- list.files(file.path(projectPath,"faustData","sampleData"))
    allPhenotypeList <- c()
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
        #
        #set up a list to store lookups
        #
        markerPhenoLookups <- vector(mode="list",(length(targetExprs)+1))
        names(markerPhenoLookups) <- c(names(targetExprs),"obsInSample")
        markerPhenoLookups[["obsInSample"]] <- nrow(amIn)
        #
        #get lookups for the markers that constitute the phenotype
        #
        for (mn in names(targetExprs)) {
            currentBM <- bitwhich(maxindex=nrow(amIn),x=FALSE)
            currentLookup <- which(amIn[,mn]==as.numeric(targetExprs[[mn]]))
            currentBM[currentLookup] <- TRUE
            markerPhenoLookups[[mn]] <- currentBM
        }
        #
        #record data from sample
        #
        allPhenotypeList <- append(allPhenotypeList,list(markerPhenoLookups))
        names(allPhenotypeList)[length(allPhenotypeList)] <- sampleName
    }
    return(allPhenotypeList)
}

.getPhenotypeCountDF <- function(projectPath,
                                markersInPheno,
                                targetExprs,
                                selectedChannels,
                                phenotypeLookupList,
                                markersInParentPhenotype)
{
    allSampleNames <- list.files(file.path(projectPath,"faustData","sampleData"))
    pcVec <- ccVec <- snVec <- c()
    for (sampleName in allSampleNames) {
        #
        #get counts for the markers that constitute the phenotype
        #
        obsInSample <- phenotypeLookupList[[sampleName]][["obsInSample"]]
        childBM <- bitwhich(maxindex=obsInSample,x=TRUE)
        for (mn in markersInPheno) {
            currentBM <- phenotypeLookupList[[sampleName]][[mn]]
            childBM <- (childBM & currentBM) #S3 method from bit library
        }
        childCount <- length(which(as.logical(childBM)))
        if (is.null(markersInParentPhenotype)) {
            parentCount <- obsInSample
        }
        else {
            parentBM <- bitwhich(maxindex=obsInSample,x=TRUE)
            for (mn in markersInParentPhenotype) {
                currentBM <- phenotypeLookupList[[sampleName]][[mn]]
                parentBM <- (parentBM & currentBM) #S3 method from bit library
            }
            parentCount <- length(which(as.logical(parentBM)))
        }
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

.getPvalueForMarkers <- function(projectPath,
                                candidateMarkers,
                                targetExprs,
                                selectedChannels,
                                metaDataDF,
                                modelStr,
                                phenotypeLookupList,
                                markersInParentPhenotype
                                )
{
    #
    #get counts for the phenotype consisting of the candidate markers
    #
    phenoDF <- .getPhenotypeCountDF(
        projectPath=projectPath,
        markersInPheno=candidateMarkers,
        targetExprs=targetExprs,
        selectedChannels=selectedChannels,
        phenotypeLookupList=phenotypeLookupList,
        markersInParentPhenotype=markersInParentPhenotype
    )
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
