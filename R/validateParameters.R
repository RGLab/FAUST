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
                                annotationsApproved,
                                archDescriptionList)
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
    if (!max(Reduce(union,list(a=(channelBounds==""),b=(is.matrix(channelBounds)),c=is.list(channelBounds))))) {
        print("channelBounds parameter does not conform to one of the following.")
        print("    ")
        print('First, you can set channelBounds to the empty string "".')
        print('This will cause the values to be estimated empirically.')
        print("    ")
        print("Second, you can set the channel bounds to be a 2 x length(activeChannels) matrix.")
        print("    ")
        print("Third, you can set the channel bounds to be a list of 2 x length(activeChannels) matrix.")
        print("In this third case, the length of the list must correspond to the number of unique levels")
        print("in the imputation hierarchy.")
        print("    ")
        print("Each slot in the list must be a string that maps to the level in the pData table.")
        print("The corresponding 2 x length(activeChannels) matrix will be used to determine admissible")
        print("events when densities are estimated in the annotation forest construction.")
        print("    ")
        stop("Please re-run faust with channelBounds conforming to one of these three options.")
    }
    if (is.character(channelBounds) && (channelBounds != "")) {
        print('Unsupported string passed to channelBounds parameter.')
        stop('Re-run set to channelBounds="" to compute bounds empirically.')
    }
    if (is.matrix(channelBounds)) {
        if ((nrow(channelBounds) != 2) ||
            (ncol(channelBounds) != length(activeChannels)))
        {
            stop("channelBounds must be a 2 x length(activeChannels) matrix.")
        }
        if (length(intersect(rownames(channelBounds),c("Low","High"))) != 2)
        {
            stop("row names of channelBounds must be 'Low' and 'High'.")
        }
        if (length(intersect(colnames(channelBounds),activeChannels)) != length(activeChannels))
        {
            stop("column names of channelBounds must be the channels in the activeChannels vector.")
        }
    }
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
    if ((archDescriptionList$targetArch!="slurmCluster") && (archDescriptionList$targetArch!="singleCPU")) {
        stop("Must set the targetArch slot in the parameter archDescriptionList to either 'singleCPU' or 'slurmCluster'")
    }
    if (archDescriptionList$targetArch=="slurmCluster") {
        if (length(
            setdiff(
                c("partitionID","jobPrefix","jobTime","maxNodeNum","maxTime","nodeThreadNum","sbatchFlags"),
                names(archDescriptionList))))
        {
            print("The targetArch slot of the archDescriptionList is set to 'slurmCluster'.")
            print("To use this setting, you must run faust on a system with slurm installed.")
            print("You must also provide entries for the following entries in the archDescriptionList.")
            print(c("partitionID","jobPrefix","jobTime","maxNodeNum","maxTime","nodeThreadNum","sbatchFlags"))
            stop("Please see the documentation and provide a setting for all parameters.")
        }
        if (!is.character(archDescriptionList$partitionID)) {
            stop("Please set the partitionID slot in the archDescriptionList to be a character string.")
        }
        if (!is.character(archDescriptionList$jobPrefix)) {
            stop("Please set the jobPrefix slot in the archDescriptionList to be a character string.")
        }
        if (!is.character(archDescriptionList$jobTime)) {
            stop("Please set the jobTime slot in the archDescriptionList to be a character string 'HH:MM:SS'.")
        }
        if (!is.numeric(archDescriptionList$maxNodeNum)) {
            stop("Please set the maxNodeNum slot in the archDescriptionList to be a non-negative numeric value.")
        }
        if (!is.numeric(archDescriptionList$maxTime)) {
            stop("Please set the maxTime slot in the archDescriptionList to be a non-negative numeric value.")
        }
        if (!is.character(archDescriptionList$nodeThreadNum)) {
            stop("Please set the nodeThreadNum slot in the archDescriptionList to be a character string.")
        }
        if (!is.character(archDescriptionList$sbatchFlags)) {
            stop("Please set the sbatchFlags slot in the archDescriptionList to be a character string.")
        }
    }
    return()
}

    
