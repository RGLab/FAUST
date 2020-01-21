.initializeFaustDataDir <- function(
                                    projectPath,
                                    activeChannels,
                                    channelBounds,
                                    startingCellPop
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
        saveRDS(startingCellPop,
                file.path(normalizePath(projectPath),
                          "faustData",
                          "metaData",
                          "channelBounds.rds"))
    }
    return()
}
