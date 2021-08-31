#' Faster Annotation Using Shape-constrained Trees
#'
#' growAnnotationForest is used to produce an annotation forest for a dataSet.
#' 
#' @param dataSet The input dataSet to be clustered by faust.
#'   Must be an R matrix. The dataSet should only contain numerical data.
#' 
#' @param numberIterations The number of faust iterations to perform.
#' 
#' @param pValueThreshold The significance level for the DIP test of Hartigan
#'   and Hartigan (1985). faust uses this level to choose to partition a projection
#'   in the dataSet with the taut-string of Davies and Kovac (2004). Setting it
#'   to smaller values (such as a conventional 0.05) usually produces faster
#'   exaustive clusterings of the dataSet.
#'
#' @param minimumClusterSize The size of the smallest admissible candidate cluster.
#'   The recursive search for candidate clusters will not consider subsets of observations
#'   below this value as admissible candidates -- it terminates along such branches without
#'   recording the subsets encountered. Increasing this parameter usually increases
#'   the speed of a single faust iteration.
#'
#' @param maximumSearchDepth An upper bound on the depth of a search tree for candidate clusters.
#'   The recursive search terminates if it exceeds this depth. Since the annotation forest
#'   is contructing boundaries relative to a *population*, it is set to a low default number.
#' 
#'   This is done under the assumption that a search path that terminates far down a gating tree
#'   is potentially finding differences that are not reflective of population boundaries. 
#'
#' @param maximumGatingNum An upper bound on the number of gating examples in the annotatino forest.
#'   If a faust iteration is exhaustively searching the space of clustering trees, the search will
#'   terminate if faust finds a collection of candidate clusters that exceed this bound.
#'
#'   If a faust iteration performs a randomCandidateSearch (or if any values are restricted in the search
#'   process using the resValMatrix) \emph{\bold{the faust search will continue until it
#'   finds this number of gating examples}}.
#'
#' @param maximumAntimodeNumber An upper bound on the number of possible antimodes found by the
#'   taut-string estimate of a projection. If the taut-string produces more than the
#'   maximumAntimodeNumber number of antimodes, the recursive search terminates along that branch.
#'
#' @param maximumSearchTime The maximum time in seconds a faust iteration can search for candidate
#'   clusters. If a single iteration exceeds this value in either the intial search phase or
#'   and residual search phases, that faust iteration is abandonded and the next iteration is 
#'   attempted.
#'
#'   On large dataSets (with either a large number of observations or variables), the
#'   initial addition of noise can produce configurations of the dataSet which take much longer
#'   to cluster than the average faust iteration. This parameter is useful to prevent such iterations
#'   from extending a single faust call from running too long.
#'
#' @param maximumRunTime The total time in seconds allowed for the growth of the annotation forest.
#' If this time is exceeding during a FAUST iteration, the iteration will take as long as needed to finish
#' before the annotation forest is returned to the user.
#' 
#' @param annotationVec A vector of strings that will be used to annotate the components of the annotation forest.
#'   For the faust clustering to be interpetable, meaningful strings must be provided here.
#'
#' @param numberOfThreads The addition of noise to the dataSet and the random search for gate locations
#'   are parallelized. For numberOfThreads greater than 2, numberOfThreads-1 worker threads will add
#'   noise to the dataSet, conduct the search for candidate clusters, and annotate the selected clusters.
#'   The default value of 0 indicates faust should detect the maximum number of threads supported by the hardware
#'   and then use that value.
#' 
#' @param randomCandidateSearch Boolean flag. If true, the initial search for candidate clusters will
#'   sample from the space of clustering trees. At each node in the clustering tree, faust will pick 
#'   the next clustering branch uniformly at random the set of admissible projections (according to the dip).
#'   Especially in small dataSets, this can cause the faust search to find the same clustering tree multiple
#'   times.
#'
#' @param anyValueRestricted Boolean flag. If true, indicates the user has some prior knowledge about bounds of
#'  admissible values. Such values must be identified by coordinate projection in the resValMatrix.
#'
#'  The motivation for this parameter is that some technologies produce dataSets which are: very zero-inflated;
#'  have a minimum limit of detection; have an upper limit of detection. Hence, one might wish to indicate that
#'  all values below some threshold be treated as low \emph{by default}, or above a threshold be treated as high
#'  \emph{by default}. Of critical importance: restricted values \emph{\bold{are not used by the taut-string density estimator
#'  in the search for candidate clusters}}.
#'  
#' @param resValMatrix R matrix indicating restricted values. Contains only three entries: {0,1,2}. If entry (i,j) is 0, it indicates
#'  that observation i is unrestricted along coordinate j (and so is used by the taut-string in the search for cnadidate clusters).
#'
#'  If entry (i,j) is 1, it indicates observation i is \emph{restricted from below} along coordinate j: this observation will 
#'  automatically be treated as "Lowest" along coordinate j, and will not be used by the taut-string in the search for candidate clusters.
#'
#'  If entry (i,j) is 2, it indicates observation i is \emph{restricted from above} along coordinate j: this observation will 
#'  automatically be treated as "Highest" along coordinate j, and will not be used by the taut-string in the search for candidate clusters.
#'
#' @param getDebugInfo Boolean flag. If true, debugging info pertaining to the annotation forest is printed to file.
#'
#' @param cutPointUpperBound Integer upper bound on the number of cut points that will be collected
#' for a column projection during the construction of the annotation forest. If, during the search for candidate clusters,
#' faust partitions a column projection into a number of modal groups that exceeds this number, the search will proceed
#' but the location of the cut points will not be added to the annotation forest. 
#' 
#' @param gaussianScaleParameter Noise is added in a two-step procedure for FAUST. First, uniform noise is added to
#' break ties observed in the columns of the dataSet. After uniform noise is added, Gaussian noise is added to perturb
#' the relative order of the observations in each column of the dataSet. The standard-devaition of the Gaussian noise
#' depends on the distance to the neighboring order statistics of an observation.
#'
#' The default value of 4.0 preserves the relative ordering of observations in a column, with high probability. As the
#' value decreases, observations are more and more likely to switch position in the order-statistics. 
#'
#' @param randomSeed Set to a large integer value to reproduce faust runs when numberOfThreads==1.
#' Default value of 0 leads to non-reproducible FAUST runs regardless of numberOfThreads value.
#'
#' @param allowRepeatedSplitting Set to FALSE by default. If TRUE, indicates the search for candidate clusters can
#' induce modal clusters along the same column vector of dataSet multiple times throughout the search. No matter
#' the value of this parameter (TRUE OR FALSE), the search for candidate clusters is never allowed to induce modal clusters
#' along the same column vector twice in a row: from the point of view of FAUST, this can only arise due to technical
#' error.
#' 
#' @param subSamplingThreshold Determines how large a data matrix needs to be in order to start sub-sampling for density
#' estimation.
#'
#' @param subSampleSize Number of observations a conditional density will sample.
#'
#' @param subSampleIter Number of sub-sampling iterations.
#' 
#' @param recordCounts Record counts of which observations are split at different levels of annotation forest.
#' Memory intensive. In runs using multiple threads, each thread maintains a local copy of a count matrix with
#' the same dimensions as the data set, and counts are copied to the master thread. 
#' 
#' @param recordIndices Record boolean vector of indices.
#' Memory intensive. Used for conditional plotting.
#' 
#' @return growAnnotationForest returns an R list with ncol(dataSet) entires. Each entry corresponds to a column in dataSet.
#' Each entry in the list is a vector of cut points placed by faust along that column projection in the search for
#' candidate clusters. The distribution of cut points can be used to set annotation boundaries for the
#' annForestVals in the faust function.
#'
#' @examples
#' clusterMatrix <- as.matrix(iris[,-5])
#' annotationForest <- growAnnotationForest(dataSet=clusterMatrix,
#'                          numberIterations=1000)
#' @export
#' @md

growAnnotationForest <- function(dataSet,
                                 numberIterations = 1,
                                 pValueThreshold = 0.35,
                                 minimumClusterSize = 25,
                                 maximumSearchDepth = ncol(dataSet),
                                 maximumGatingNum = 1000,
                                 maximumAntimodeNumber = 100,
                                 maximumSearchTime=1e6,
                                 maximumRunTime = 1e6,
                                 annotationVec = colnames(dataSet),
                                 numberOfThreads=1,
                                 randomCandidateSearch=FALSE,
                                 anyValueRestricted=FALSE,
                                 resValMatrix=matrix(0,nrow=2,ncol=2),
                                 getDebugInfo=FALSE,
                                 cutPointUpperBound=2,
                                 gaussianScaleParameter=4,
                                 randomSeed=0,
                                 allowRepeatedSplitting=FALSE,
                                 subSamplingThreshold=1e6,
                                 subSampleSize=1e6,
                                 subSampleIter=1,
                                 recordCounts=FALSE,
                                 recordIndices=FALSE){
    if (!is.matrix(dataSet))
        stop("Must provide a numeric R matrix")
    if (!is.double(pValueThreshold))
        stop("The p-value threshold must be a double precison number.")
    if ((pValueThreshold < 0.01) || (pValueThreshold > 0.99))
        stop("The p-value threshold for the dip test must be between  0.01 and 0.99")
    if ((!is.numeric(minimumClusterSize)) || (minimumClusterSize < 10))
        stop("A cluster cannot be smaller than 10 observations")
    if ((!is.numeric(maximumSearchDepth)) || (maximumSearchDepth < 1))
        stop("Search depth must be a natural number > 1")
    if ((!is.numeric(maximumGatingNum)) || (maximumGatingNum < 2))
        stop("Must search for more than 1 gating example")
    if (nrow(dataSet) < minimumClusterSize)
        stop("Too few observations for data")
    if (max(is.na(annotationVec)))
        stop('User must assign column names to input data set. colnames(dataSet) <- c("user", "entry",...)')
    if (maximumRunTime <= 0)
        stop("User must assign positive value for maximumRunTime (in seconds)")
    if (gaussianScaleParameter < 0.1)
        stop("Gaussian scale parameter smaller than 0. Not possible.")
    if (!(identical(randomCandidateSearch,TRUE) || identical(randomCandidateSearch,FALSE))) 
        stop("User must set randomCandidateSearch either to TRUE or FALSE.")
    if (!(identical(anyValueRestricted,TRUE) || identical(anyValueRestricted,FALSE))) 
        stop("User must set anyValueRestricted either to TRUE or FALSE.")
    if (!(identical(allowRepeatedSplitting,TRUE) || identical(allowRepeatedSplitting,FALSE))) 
        stop("User must set allowRepeatedSplitting either to TRUE or FALSE.")


    aCols <- apply(dataSet,2,length)
    uCols <- apply(dataSet,2,function(x){length(unique(x))})
    rCols <- uCols/aCols
    mCol <- min(rCols)

    firstIteration <- TRUE
    remainingIterations <- numberIterations
    startTime <- proc.time()
    timeDiff <- proc.time()-startTime
    elapsedTime <- as.numeric(timeDiff[3])
    if (randomSeed > 0) {
        set.seed(randomSeed)
    }
    while ((remainingIterations > 0) && (elapsedTime <= maximumRunTime)) {
        if (((remainingIterations %% 100) == 0) || (getDebugInfo)) {
            print(paste("remaining iterations: ", remainingIterations, sep=""))
        }
        if (randomSeed > 0) {
            seedVal <- ceiling(stats::runif(1,min=10,max=100000000))
        }
        else {
            seedVal <- 0
        }
        gatingLocs <- cppGrowAnnotationForest(dataSet,
                                              as.double(pValueThreshold),
                                              as.integer(minimumClusterSize),
                                              allowRepeatedSplitting,
                                              maximumSearchDepth,
                                              maximumGatingNum,
                                              getDebugInfo,
                                              maximumAntimodeNumber,
                                              randomCandidateSearch,
                                              numberOfThreads,
                                              anyValueRestricted,
                                              resValMatrix,
                                              as.integer((cutPointUpperBound+2)), #add two for c++ representation
                                              maximumSearchTime,
                                              gaussianScaleParameter,
                                              seedVal,
                                              subSamplingThreshold,
                                              subSampleSize,
                                              subSampleIter,
                                              recordCounts,
                                              recordIndices)
        if (firstIteration) {
            gatingForest <- gatingLocs
            firstIteration <- FALSE
        }
        else {
            #append results to iterated runs.
            for (i in seq(length(gatingForest[["gateData"]]))) {
                gatingForest[["gateData"]][[i]] <- append(gatingForest[["gateData"]][[i]],
                                                          gatingLocs[["gateData"]][[i]])
            }
            gatingForest[["subsetCounts"]] <-  gatingForest[["subsetCounts"]] + gatingLocs[["subsetCounts"]]
            gatingForest[["subsetDenom"]] <-  gatingForest[["subsetDenom"]] + gatingLocs[["subsetDenom"]]
            for (i in seq(length(gatingForest[["indexData"]]))) {
                gatingForest[["indexData"]][[i]] <- append(gatingForest[["indexData"]][[i]],
                                                           gatingLocs[["indexData"]][[i]])
            }
            for (i in seq(length(gatingForest[["indexDepthData"]]))) {
                gatingForest[["indexDepthData"]][[i]] <- append(gatingForest[["indexDepthData"]][[i]],
                                                                gatingLocs[["indexDepthData"]][[i]])
            }
        }
        #decrement the number of iterations and compute time elapsed.
        remainingIterations <- (remainingIterations - 1)
        timeDiff <- proc.time()-startTime
        elapsedTime <- as.numeric(timeDiff[3])
    }
    #print(paste("Completed ", (numberIterations-remainingIterations), " iterations.",sep=""))
    afListEntries <- c("cutPoints","numCuts","cutDepth","nodePathScore","nodePopSize")
    finalAnnotation <- c()
    for (anColNum in seq(length(annotationVec))) {
        descUpdate <- paste0(annotationVec[anColNum],"_",afListEntries)
        finalAnnotation <- append(finalAnnotation,descUpdate)
    }
    names(gatingForest[["gateData"]]) <- finalAnnotation
    colnames(gatingForest[["subsetCounts"]]) <- annotationVec
    names(gatingForest[["subsetDenom"]]) <- annotationVec
    names(gatingForest[["indexData"]]) <- annotationVec
    names(gatingForest[["indexDepthData"]]) <- annotationVec
    return(gatingForest)
}
