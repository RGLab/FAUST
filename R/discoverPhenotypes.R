discoverPhenotypes <- function(projectPath,
                               numScampIter,
                               nameOccuranceNum,
                               debugFlag,
                               threadNum,
                               seedValue,
                               archDescriptionList,
                               plottingDevice)
{
    #test parameters for validity. stop faust run if invalid settings detected.
    .validateDiscoveryParameters(
        projectPath = projectPath,
        debugFlag = debugFlag,
        threadNum = threadNum,
        seedValue = seedValue,
        archDescriptionList=archDescriptionList
    )

    if (debugFlag) print("Discovering phenotypes across experimental units.")
    .clusterLevelsWithScamp(
        projectPath = projectPath,
        numScampIter = numScampIter,
        nameOccuranceNum = nameOccuranceNum,
        debugFlag = debugFlag,
        threadNum = threadNum,
        seedValue = seedValue,
        archDescriptionList = archDescriptionList
    )

    .plotPhenotypeFilter(
        projectPath=projectPath,
        nameOccuranceNum=nameOccuranceNum,
        plottingDevice=plottingDevice
    )

    if (debugFlag) print("Gating populations.")
    .gateScampClusters(
        projectPath = projectPath,
        debugFlag = debugFlag
    )

    if (debugFlag) print("Generating faust count matrix.")
    .getFaustCountMatrix(
        projectPath = projectPath,
        debugFlag = debugFlag
    )
}
