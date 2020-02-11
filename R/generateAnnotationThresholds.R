generateAnnotationThresholds <- function(gatingSet,
                                         projectPath,
                                         experimentalUnit,
                                         imputationHierarchy,
                                         activeChannels,
                                         channelBounds,
                                         startingCellPop,
                                         numForestIter,
                                         depthScoreThreshold,
                                         selectionQuantile,
                                         seedValue,
                                         threadNum,
                                         debugFlag,
                                         supervisedList,
                                         archDescriptionList,
                                         annotationsApproved,
                                         drawAnnotationHistograms,
                                         plottingDevice)
{
    #first, test parameters for validity. stop faust run if invalid settings detected.
    .validateForestParameters(
        activeChannels = activeChannels,
        channelBounds = channelBounds,
        startingCellPop = startingCellPop,
        projectPath = projectPath,
        depthScoreThreshold = depthScoreThreshold,
        selectionQuantile = selectionQuantile,
        debugFlag = debugFlag,
        threadNum = threadNum,
        seedValue = seedValue,
        supervisedList = supervisedList,
        annotationsApproved = annotationsApproved,
        archDescriptionList=archDescriptionList
    )


    #set up the faustData directory for check-pointing/metadata storage.
    .initializeFaustDataDir(
        projectPath = projectPath,
        activeChannels = activeChannels,
        channelBounds = channelBounds,
        startingCellPop = startingCellPop,
        supervisedList=supervisedList
    )

    #construct the analysis map using the metadata stored in the gating set.
    #the analysis map is a data frame that links samples to their experimental unit
    #and their level in the imputation hierarchy.
    .constructAnalysisMap(
        projectPath = projectPath,
        gspData = flowWorkspace::pData(gatingSet),
        sampNames = flowWorkspace::sampleNames(gatingSet),
        experimentalUnit = experimentalUnit,
        imputationHierarchy = imputationHierarchy,
        debugFlag = debugFlag
    )

    #begin method processing. copy data to projectPath from gatingSet.
    if (debugFlag) print("Begin data extraction.")
    .extractDataFromGS(
        gs = gatingSet,
        activeChannels = activeChannels,
        startingCellPop = startingCellPop,
        projectPath = projectPath,
        debugFlag = debugFlag
    )

    #make sure the channel bounds conform to internal requirements
    #test to see if there have been changes between faust runs.
    .processChannelBounds(
        samplesInExp = flowWorkspace::sampleNames(gatingSet),
        projectPath = projectPath,
        channelBounds = channelBounds,
        debugFlag = debugFlag
    )

    if (debugFlag) print("Making restriction matrices.")
    .makeRestrictionMatrices(
        samplesInExp = flowWorkspace::sampleNames(gatingSet),
        channelBounds = channelBounds,
        projectPath = projectPath,
        debugFlag = debugFlag
    )

    #accumulate data into the experimental units.
    if (debugFlag) print("Collecting data into experimental units.")
    .prepareFirstAL(projectPath = projectPath)


    #
    #ALL uses of startingCellPop must be updated to sanitizedCellPopStr
    #
    #start the annotation process
    if (!file.exists(file.path(normalizePath(projectPath),
                               "faustData",
                               "metaData",
                               "bigForestDone.rds"))) {
        #in large experiments, this can be a costly step without sub-sampling.
        #often we will want to supervise the results after growing the forest,
        #so only grow it on an as-need basis
        .growAnnForest(
            activeChannels = activeChannels,
            numIter = numForestIter,
            debugFlag = debugFlag,
            threadNum = threadNum,
            seedValue = seedValue,
            projectPath = projectPath,
            archDescriptionList = archDescriptionList
        )
        bigForestDone <- TRUE
        saveRDS(bigForestDone,
                file.path(normalizePath(projectPath),
                          "faustData",
                          "metaData",
                          "bigForestDone.rds"))
    }

    if (debugFlag) print("Selecting standard set of channels across experiment using depth score.")
    .selectChannels(
        depthScoreThreshold = depthScoreThreshold,
        selectionQuantile = selectionQuantile,
        projectPath = projectPath
    )
    #selC no longer returned
    #must map to reads from metaData/initSelC.rds

    if (debugFlag) print("Reconciling annotation boundaries across experiment.")
    .reconcileAnnotationBoundaries(
        projectPath = projectPath,
        debugFlag = debugFlag
    )

    .superviseReconciliation(
        projectPath = projectPath,
        debugFlag = debugFlag
    )


    if (debugFlag) print("Writing annotation matrices to file.")
    .mkAnnMats(
        projectPath = projectPath
    )

    if (debugFlag) print("Generating depth score plot.")
    .plotScoreLines(
        projectPath = projectPath,
        depthScoreThreshold = depthScoreThreshold,
        selectionQuantile = selectionQuantile,
        plottingDevice = plottingDevice
    )

    if (debugFlag) print("Generating marker boundary histograms.")
    .plotMarkerHistograms(
        projectPath = projectPath,
        plottingDevice = plottingDevice
    )

    if (drawAnnotationHistograms) {
        if (debugFlag) print("Generating annotation boundary histograms.")
        .plotSampleHistograms(
            projectPath = projectPath,
            plottingDevice=plottingDevice
        )
    }

    print("********************************************************")
    print("FAUST has selected a subset of the marker panel using")
    print("the depth score at the specified selection quantile.")
    print("Annotation thresholds have been generated for all selected markers.")
    print(paste0("Plots have been written to file in the directory ",
                 paste0(file.path(normalizePath(projectPath),
                                  "faustData",
                                  "plotData"))))
    print("")
    print("Review these plots to ensure all desired markers have been selected.")
    print("If too many/too few markers have been selected, modify the parameters")
    print("depthScoreThreshold and selectionQuantile.")
    print("")
    print(paste0("Also review the annotation thresholds displayed on the sample-level histograms in",
                 paste0(file.path(normalizePath(projectPath),
                                  "faustData",
                                  "plotData"))))
    print("If you wish to modify the placement of the annotation thresholds,")
    print("change the parameters supervisedList.")
    print("")
    print("Changing the Low/High values in the channelBounds matrix will also affect")
    print("the placement of annotation thresholds. It is the most effective way to directly modify")
    print("their placement. However, when you modify the Low/High values in the channelBounds matrix,")
    print("the FAUST method will regrow the entire annotation forest.")
    print("")
    print("Once you are satisfied with the annotation boundary")
    print("placement, set the parameter")
    print("annotationsApproved=TRUE in the faust R function.")
    print("FAUST will then conduct phenotype discovery and ")
    print("annotation on each experimental unit.")
    print("********************************************************")
    if (annotationsApproved)
    {
        saveRDS(annotationsApproved,
                file.path(normalizePath(projectPath),
                          "faustData",
                          "metaData",
                          "annotationsApproved.rds"))

    }
    return()
}
