.initializeFaustDataDir <- function(
                                    projectPath,
                                    activeChannels,
                                    channelBounds,
                                    startingCellPop,
                                    supervisedList
                                    )
{
    #initialize pipeline directories. save metadata.
    if (!dir.exists(file.path(normalizePath(projectPath),
                              "faustData")))
    {
        dir.create(file.path(normalizePath(projectPath),
                             "faustData"))
    }
    if (!dir.exists(file.path(normalizePath(projectPath),
                              "faustData",
                              "sampleData")))
    {
        dir.create(file.path(normalizePath(projectPath),
                             "faustData",
                             "sampleData"))
    }
    if (!dir.exists(file.path(normalizePath(projectPath),
                              "faustData",
                              "metaData")))
    {
        dir.create(file.path(normalizePath(projectPath),
                             "faustData",
                             "metaData"))
    }
    if (!dir.exists(file.path(normalizePath(projectPath),
                              "faustData",
                              "plotData")))
    {
        dir.create(file.path(normalizePath(projectPath),
                             "faustData",
                             "plotData"))
    }
    if (!dir.exists(file.path(normalizePath(projectPath),
                              "faustData",
                              "plotData",
                              "histograms")))
    {
        dir.create(file.path(normalizePath(projectPath),
                             "faustData",
                             "plotData",
                             "histograms"))
    }
    if (!file.exists(file.path(normalizePath(projectPath),
                               "faustData",
                               "metaData",
                               "activeChannels.rds")))
    {
        saveRDS(activeChannels,
                file.path(normalizePath(projectPath),
                          "faustData",
                          "metaData",
                          "activeChannels.rds"))
    }
    if (!file.exists(file.path(normalizePath(projectPath),
                               "faustData",
                               "metaData",
                               "startingCellPop.rds")))
    {
        saveRDS(startingCellPop,
                file.path(normalizePath(projectPath),
                          "faustData",
                          "metaData",
                          "startingCellPop.rds"))
    }
    if (!file.exists(file.path(normalizePath(projectPath),
                               "faustData",
                               "metaData",
                               "channelBounds.rds")))
    {
        saveRDS(channelBounds,
                file.path(normalizePath(projectPath),
                          "faustData",
                          "metaData",
                          "channelBounds.rds"))
    }

    #always sanitize the starting cell pop for problem characters.
    sanitizedCellPopStr <- gsub("[[:punct:]]","",startingCellPop)
    sanitizedCellPopStr <- gsub("[[:space:]]","",sanitizedCellPopStr)
    sanitizedCellPopStr <- gsub("[[:cntrl:]]","",sanitizedCellPopStr)
    saveRDS(sanitizedCellPopStr,
            file.path(normalizePath(projectPath),
                      "faustData",
                      "metaData",
                      "sanitizedCellPopStr.rds"))

    #always update the supervision artifacts in the metaData direcotry
    forceList <- selectionList <- preferenceList <- list()
    if (!is.na(supervisedList)) {
        #supervisedList is a named list of lists
        #name of slot in list: marker
        #list under marker slot 1: string describing type of supervision.
        #list under marker slot 2: vector of ints dictating supervision action.
        supervisedMarkers <- names(supervisedList)
        for (markerNum in seq(length(supervisedMarkers))) {
            marker <- supervisedMarkers[markerNum]
            markerList <- supervisedList[[markerNum]]
            actionType <- markerList[[1]]
            action <- markerList[[2]]
            if (actionType == "Preference")  {
                preferenceList <- append(preferenceList,list(action))
                names(preferenceList)[length(preferenceList)] <- marker
            }
            else if (actionType == "PostSelection")  {
                selectionList <- append(selectionList,list(action))
                names(selectionList)[length(selectionList)] <- marker
            }
            else if (actionType == "Force") {
                forceList <- append(forceList,list(action))
                names(forceList)[length(forceList)] <- marker
            }
            else {
                print(paste0("Requested unsupported supervision type for marker ", marker))
                print("Only 'Force', 'Preference' and 'PostSelection' supervision types are currently supported.")
                print(paste0("Ignoring requested action: ",actionType))
            }
        }
    }

    saveRDS(forceList,
            file.path(normalizePath(projectPath),
                      "faustData",
                      "metaData",
                      "forceList.rds"))

    saveRDS(selectionList,
            file.path(normalizePath(projectPath),
                      "faustData",
                      "metaData",
                      "selectionList.rds"))

    saveRDS(preferenceList,
            file.path(normalizePath(projectPath),
                      "faustData",
                      "metaData",
                      "preferenceList.rds"))

    return()
}
