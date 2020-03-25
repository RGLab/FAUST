#' Full Annotation Using Shape-constrained Trees
#'
#' Once FAUST has been used to generate a count matrix, this function
#' can be used to create a sample-by-phenotype matrix of median marker intensities
#' for a target marker of interest. The target marker must be one of the active
#' channels in the initial analysis.
#'
#' @param targetMarker The marker of interest. Function will construct a matrix
#' of median marker for this target marker, assuming it was an active channel in the
#' initial faust analyiss.
#'
#' @param projectPath The filepath locating the faustData directory of the completed
#' faust analysis.
#'
#' @param debugFlag Bolean value. Set to TRUE to print sample processing info.
#'
#' @return A sample-by-phenotype matrix of median marker intensities for the
#' target marker.
#'
#' @examples
#'
#' #getMMIForTargetMarker(targetMarker="IFNg",projectPath="/home/user",debugFlag=TRUE)
#' 
#' @export
#' @md
getMMIForTargetMarker <- function(
                                  targetMarker,
                                  projectPath=normalizePath("."),
                                  debugFlag=FALSE
                                  )
{
    initAC <- readRDS(file.path(normalizePath(projectPath),
                                "faustData",
                                "metaData",
                                "activeChannels.rds"))
    if (!(targetMarker %in% initAC)) {
        stop(paste0("targetMarker ",targetMarker," not in initial active channels."))
    }
    cnMap <- readRDS(file.path(normalizePath(projectPath),
                               "faustData",
                               "metaData",
                               "colNameMap.rds"))
    numClusters <- nrow(cnMap)
    allSampleNames <- list.files(file.path(normalizePath(projectPath),
                                           "faustData",
                                           "sampleData"))
    mmiMatrix <- matrix(NA,nrow=length(allSampleNames),ncol=numClusters)
    colnames(mmiMatrix) <- cnMap[,"newColNames",drop=TRUE]
    rownames(mmiMatrix) <- allSampleNames
    for (sn in allSampleNames) {
        if (debugFlag) print(paste0("Processing sample ",sn))
        exprsMat <- readRDS(file.path(normalizePath(projectPath),
                                      "faustData",
                                      "sampleData",
                                      sn,
                                      "exprsMat.rds"))
        fa <- read.table(file.path(normalizePath(projectPath),
                                   "faustData",
                                   "sampleData",
                                   sn,
                                   "faustAnnotation.csv"),
                         header=F,
                         sep="`",
                         stringsAsFactors=FALSE)[,1]
        for (cNum in seq(numClusters)) {
            cEnc <- cnMap[cNum,"faustColNames",drop=TRUE]
            cName <- cnMap[cNum,"newColNames",drop=TRUE]
            cLookup <- which(fa==cEnc)
            if (length(cLookup)) {
                mmi <- median(exprsMat[cLookup,targetMarker])
                mmiMatrix[sn,cName] <- mmi
            }
        }
    }
    return(mmiMatrix)
}

