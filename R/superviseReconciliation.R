.superviseReconciliation <- function(supervisedList,parentNode,projectPath=".")
{
    resListPrep <- readRDS(file.path(normalizePath(projectPath),
                                     "faustData",
                                     "gateData",
                                     paste0(parentNode,"_resListPrep.rds")))
    outList <- resListPrep
    selectedChannels <- names(resListPrep)
    supervisedChannels <- names(supervisedList)
    if (length(setdiff(supervisedChannels,selectedChannels))) {
        print("The following unselected channels (by depth score) are detected.")
        print(setdiff(supervisedChannels,selectedChannels))
        print("Proceding as if these are controlled values.")
    }
    for (channel in supervisedChannels) {
        tmpList <- outList[[channel]]
        supervision <- supervisedList[[channel]]
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
