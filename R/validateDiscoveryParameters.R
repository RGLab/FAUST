.validateDiscoveryParameters <- function(projectPath,
                                         debugFlag,
                                         threadNum,
                                         seedValue,
                                         archDescriptionList)
{
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
    if (length(projectPath) != 1)
    {
        stop("projectPath must be a single character string.")
    }
    #check to make sure parameters conform to type and content expectations.
    if (!is.logical(debugFlag))
    {
        stop("debugFlag must be set to TRUE xor FALSE.")
    }
    if ((!is.numeric(threadNum)) || (threadNum <= 0))
    {
        stop("threadNum must be an integer value larger than 0.")
    }
    if ((!is.numeric(seedValue)) || (seedValue <= 0))
    {
        stop("seedValue must be an integer value larger than 0.")
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


