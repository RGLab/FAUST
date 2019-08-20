.validateParameters <- function(activeChannels,
                                channelBounds,
                                startingCellPop,
                                projectPath,
                                depthScoreThreshold,
                                selectionQuantile,
                                debugFlag,
                                threadNum,
                                seedValue,
                                numForestIter,
                                numScampIter,
                                supervisedList,
                                annotationsApproved)
{
    #check to make sure parameters are the expected length.
    if (length(annotationsApproved) != 1)
    {
        stop("annotationsApproved must be set to a single boolean value.")
    }
    if ((!is.na(supervisedList)) && (length(setdiff(names(supervisedList),activeChannels))))
    {
        print("supervisedList is attempting to supervise a channel")
        stop("not in the activeChannels vector.")
    }
    if ((!is.na(supervisedList)) && (length(supervisedList) > length(activeChannels)))
    {
        print("supervisedList is attempting to supervise more channels")
        stop("than those listed in the activeChannels vector.")
    }
    if (length(numForestIter) != 1)
    {
        stop("numForestIter must be a single integer value greater than 1.")
    }
    if (length(numScampIter) != 1)
    {
        stop("numScampIter must be a single integer value greater than 1.")
    }
    if (length(seedValue) != 1)
    {
        stop("seedValue must be a single integer value.")
    }
    if (length(threadNum) != 1)
    {
        stop("threadNum must be a single integer value greater than 0.")
    }
    if (length(debugFlag) != 1)
    {
        stop("debugFlag must be a single boolean value.")
    }
    if (length(selectionQuantile) != 1)
    {
        stop("selectionQuantile must be a single numeric value between 0 and 1.")
    }
    if (length(depthScoreThreshold) != 1)
    {
        stop("depthScoreThreshold must be a single numeric value between 0 and 1.")
    }
    if (length(projectPath) != 1)
    {
        stop("projectPath must be a single character string.")
    }
    if (length(startingCellPop) != 1)
    {
        stop("startingCellPop must be a single character string.")
    }
        
    #check to make sure parameters conform to type and content expectations.
    if (!is.character(activeChannels)) {
        print("activeChannels must be a character vector listing the")
        stop("channels you wish you use in the pipeline.")
    }
    if (length(activeChannels) != length(unique(activeChannels))) {
        stop("activeChannels must contain each active channel in the experiment only once.")
    }
    #if (!is.matrix(channelBounds)) {
    #    stop("channelBounds must be a 2 x length(activeChannels) matrix.")
    #}
    #if ((nrow(channelBounds) != 2) ||
    #    (ncol(channelBounds) != length(activeChannels)))
    #{
    #    stop("channelBounds must be a 2 x length(activeChannels) matrix.")
    #}
    #if (length(intersect(rownames(channelBounds),c("Low","High"))) != 2)
    #{
    #    stop("row names of channelBounds must be 'Low' and 'High'.")
    #}
    #if (length(intersect(colnames(channelBounds),activeChannels)) != length(activeChannels))
    #{
    #    stop("column names of channelBounds must be the channels in the activeChannels vector.")
    #}
    if (!is.logical(annotationsApproved))
    {
        stop("annotationsApproved must be set to TRUE xor FALSE.")
    }
    if (!is.logical(debugFlag))
    {
        stop("debugFlag must be set to TRUE xor FALSE.")
    }
    if ((!is.na(supervisedList)) && (!is.list(supervisedList)))
    {
        stop("supervisedList must be a named list.")
    }
    if ((!is.numeric(threadNum)) || (threadNum <= 0))
    {
        stop("threadNum must be an integer value larger than 0.")
    }
    if ((!is.numeric(seedValue)) || (seedValue <= 0))
    {
        stop("seedValue must be an integer value larger than 0.")
    }
    if ((!is.numeric(numForestIter)) || (numForestIter < 1))
    {
        stop("numForestIter must be an integer value larger than 0.")
    }
    if ((!is.numeric(numScampIter)) || (numScampIter < 1))
    {
        stop("numScampIter must be an integer value larger than 0.")
    }
    if ((!is.numeric(depthScoreThreshold)) ||
        (depthScoreThreshold < 0) ||
        (depthScoreThreshold > 1))
    {
        stop("depthScoreThreshold must be a single numeric value between 0 and 1.")
    }
    if ((!is.numeric(selectionQuantile)) ||
        (selectionQuantile < 0) ||
        (selectionQuantile > 1))
    {
        stop("selectionQuantile must be a single numeric value between 0 and 1.")
    }

    #initialize pipeline directories. save metadata.
    if (!dir.exists(paste0(projectPath,"/faustData")))
    {
        dir.create(paste0(projectPath,"/faustData"))
    }
    if (!dir.exists(paste0(projectPath,"/faustData/sampleData")))
    {
        dir.create(paste0(projectPath,"/faustData/sampleData"))
    }
    if (!dir.exists(paste0(projectPath,"/faustData/metaData")))
    {
        dir.create(paste0(projectPath,"/faustData/metaData"))
    }
    if (!dir.exists(paste0(projectPath,"/faustData/plotData")))
    {
        dir.create(paste0(projectPath,"/faustData/plotData"))
    }
    if (!dir.exists(paste0(projectPath,"/faustData/plotData/histograms")))
    {
        dir.create(paste0(projectPath,"/faustData/plotData/histograms"))
    }
    if (!file.exists(paste0(projectPath,"/faustData/metaData/activeChannels.rds")))
    {
        saveRDS(activeChannels,paste0(projectPath,"/faustData/metaData/activeChannels.rds"))
    }
    if (!file.exists(paste0(projectPath,"/faustData/metaData/startingCellPop.rds")))
    {
        saveRDS(startingCellPop,paste0(projectPath,"/faustData/metaData/startingCellPop.rds"))
    }
    return()
}

    
