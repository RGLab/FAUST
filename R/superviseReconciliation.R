.superviseReconciliation <- function(projectPath,debugFlag)
{

    parentNode <- readRDS(file.path(normalizePath(projectPath),
                                    "faustData",
                                    "metaData",
                                    "sanitizedCellPopStr.rds"))
    selectionList <- readRDS(file.path(normalizePath(projectPath),
                                       "faustData",
                                       "metaData",
                                       "selectionList.rds"))
    if (length(selectionList) > 0){
        if (debugFlag) print("Selection specific reconciled annotation boundaries.")
        resListPrep <- readRDS(file.path(normalizePath(projectPath),
                                         "faustData",
                                         "gateData",
                                         paste0(parentNode,"_resListPrep.rds")))
        outList <- resListPrep
        selectedChannels <- names(resListPrep)
        supervisedChannels <- names(selectionList)
        if (length(setdiff(supervisedChannels,selectedChannels))) {
            print("The following unselected channels (by depth score) are detected.")
            print(setdiff(supervisedChannels,selectedChannels))
            print("Proceding as if these are controlled values.")
        }
        for (channel in supervisedChannels) {
            #use the selection list to set the standard
            tmpList <- outList[[channel]]
            supervision <- selectionList[[channel]]
            if (length(supervision) > length(tmpList[[1]])) {
                stop("Attempting to set more gates than exist post-reconciliation.")
            }
            if (max(supervision) > length(tmpList[[1]])) {
                stop("Attempting to set a gate beyond the last gate existing post-reconciliation.")
            }
            if (min(supervision) < 1) {
                stop("Attempting to set a gate beneath the first gate existing post-reconciliation.")
            }
            for (gateNum in seq(length(tmpList))) {
                tmpList[[gateNum]] <- tmpList[[gateNum]][supervision]
            }
            outList[[channel]] <- tmpList
        }
        saveRDS(outList,
                file.path(normalizePath(projectPath),
                          "faustData",
                          "gateData",
                          paste0(parentNode,"_resList.rds")))
    }
    else {
        file.copy(
            from = file.path(normalizePath(projectPath),
                             "faustData",
                             "gateData",
                             paste0(parentNode,"_resListPrep.rds")),
            to = file.path(normalizePath(projectPath),
                           "faustData",
                           "gateData",
                           paste0(parentNode,"_resList.rds")),
            overwrite = TRUE
        )
    }
    return()
}


