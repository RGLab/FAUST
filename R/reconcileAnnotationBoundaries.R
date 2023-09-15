.reconcileAnnotationBoundaries <- function(projectPath,
                                           debugFlag)
{
    #this function contains the logic to standardize the annotation thresholds 
    #across experimental units.
    if (!dir.exists(file.path(normalizePath(projectPath),
                              "faustData",
                              "gateData")))
    {
        dir.create(file.path(normalizePath(projectPath),
                             "faustData",
                             "gateData"))
    }

    parentNode <- readRDS(file.path(normalizePath(projectPath),
                                    "faustData",
                                    "metaData",
                                    "sanitizedCellPopStr.rds"))

    analysisMap <- readRDS(file.path(normalizePath(projectPath),
                                     "faustData",
                                     "metaData",
                                     "analysisMap.rds"))

    selectedChannels <- readRDS(file.path(normalizePath(projectPath),
                                          "faustData",
                                          "metaData",
                                          "initSelC.rds"))

    forceList <- readRDS(file.path(normalizePath(projectPath),
                                   "faustData",
                                   "metaData",
                                   "forceList.rds"))

    preferenceList <- readRDS(file.path(normalizePath(projectPath),
                                       "faustData",
                                       "metaData",
                                       "preferenceList.rds"))

    uniqueExpUnits <- unique(analysisMap[,"experimentalUnit",drop=TRUE])
    uniqueIH <- unique(analysisMap[,"impH",drop=TRUE])
    gateList <- .makeGateList(
        uniqueExpUnits=uniqueExpUnits,
        projectPath=projectPath,
        parentNode=parentNode,
        selectedChannels=selectedChannels
    )
    if (length(gateList)) {
        #the gateList is stratifed by the experimental unit.
        saveRDS(gateList,
                file.path(normalizePath(projectPath),
                          "faustData",
                          "gateData",
                          paste0(parentNode,"_rawGateList.rds")))
        #hence the number of rows in the numGateMatrix is the number of experimental units.
        numGateMatrix <- Reduce(rbind,
                                lapply(gateList,
                                       function(x){unlist(lapply(x,
                                                                 function(y){ifelse(is.na(y[1]),0,length(y))}))}))
        if (!is.matrix(numGateMatrix)) numGateMatrix <- t(as.matrix(numGateMatrix))
        rownames(numGateMatrix) <- names(gateList)
        #numSel is a numeric vector.
        #entry is the number of gates selected for the associated marker that annotates the slot.
        numSel <- .getConsistentNumberOfGatesForMarkers(
            numGateMatrix=numGateMatrix,
            preferenceList=preferenceList,
            projectPath=projectPath,
            debugFlag=debugFlag
        )
        #resList is a container for the annotation boundaries for all the markers.
        resList <- rep(list(NA),ncol(numGateMatrix))
        names(resList) <- colnames(numGateMatrix)
        for (channel in names(numSel)) {
            #for each marker, get a standard set of annotation boundaries.
            resListUpdate <- rep(list(NA),nrow(numGateMatrix))
            names(resListUpdate) <- rownames(numGateMatrix)
            #gateNumber stores the standard number of thresholds for the marker across the experiment.
            gateNumber <- as.numeric(numSel[channel])
            numLookup <- which(numGateMatrix[,channel]==gateNumber)
            matchingExpUnits <- rownames(numGateMatrix)[numLookup]
            for (currentIH in uniqueIH) {
                #for imputation hierarchy, standardize number of thresholds
                allExpUnitsInIH <- sort(unique(analysisMap[which(analysisMap[,"impH"]==currentIH),"experimentalUnit",drop=TRUE]))
                matchingExpUnitsInIH <- intersect(matchingExpUnits,allExpUnitsInIH)
                if (length(matchingExpUnitsInIH)) {
                    #there are experimental units in the imputation hierarchy that have the gateNumber of thresholds.
                    #standardize other units in the ih to these boundaries.
                    gateMatrix <- matrix(nrow=0,ncol=gateNumber)
                    for (expUnit in matchingExpUnitsInIH) {
                        gateData <- sort(gateList[[expUnit]][[channel]])
                        gateMatrix <- rbind(gateMatrix,gateData)
                        rownames(gateMatrix)[nrow(gateMatrix)] <- expUnit
                    }
                    cbLookup <- which(rownames(numGateMatrix) %in% allExpUnitsInIH)
                    possibleGates <- as.numeric(names(table(numGateMatrix[cbLookup,channel])))
                    possibleGates <- setdiff(possibleGates,c(0,gateNumber))
                    #update the resList with the thresholds
                    for (expUnit in rownames(gateMatrix)) {
                        resListUpdate[[expUnit]] <- sort(gateMatrix[which(rownames(gateMatrix)==expUnit),])
                    }
                    #check for outliers
                    mgMed <- apply(gateMatrix,2,stats::median)
                    mgMAD <- apply(gateMatrix,2,stats::mad)
                    #if we set the gate using only one example (via supervision, or due to sparsity in the
                    #level of the imputation hierarchy), the MAD is 0.
                    #set to (-Inf,Inf) because want to keep what's found by annotation forest.
                    if (any(mgMAD == 0)) {
                        mgMAD <- Inf
                    }
                    lowVal <- mgMed - (stats::qt(0.975,df=1) * mgMAD)
                    highVal <- mgMed + (stats::qt(0.975,df=1) * mgMAD)
                    chkMatrix <- Reduce(cbind,
                                        lapply(seq(ncol(gateMatrix)),
                                               function(x){(gateMatrix[,x] >= lowVal[x])+(gateMatrix[,x] <= highVal[x])}))
                    #map outliers to medians
                    amendedNames <- c()
                    for (i in seq(ncol(chkMatrix))) {
                        badLookup <- which(chkMatrix[,i] != 2)
                        if (length(badLookup)) {
                            amendedNames <- append(amendedNames,(rownames(gateMatrix)[badLookup]))
                            goodVals <- gateMatrix[-badLookup,i]
                            newVal <- stats::median(goodVals)
                            gateMatrix[badLookup,i] <- newVal
                        }
                    }
                    amendedNames <- sort(unique(amendedNames))
                    if (length(amendedNames)) {
                        for (changeName in amendedNames) {
                            resListUpdate[[changeName]] <- sort(gateMatrix[which(rownames(gateMatrix)==changeName),])
                        }
                    }
                    finalVals <- apply(gateMatrix,2,stats::median)
                    #for experimental units with different numbers of gates than the selection,
                    #impute or delete boundaries across the imputation hierarchy
                    #so they adhere to the standard number
                    if (length(possibleGates)) {
                        for (gateNum in possibleGates) {
                            modLookup <- which(numGateMatrix[,channel]==gateNum)
                            allModSamples <- rownames(numGateMatrix)[modLookup]
                            modSamples <- intersect(allModSamples,allExpUnitsInIH)
                            for (modName in modSamples) {
                                modVals <- gateList[[modName]][[channel]]
                                if (gateNum < length(finalVals)) newModVals <- .upConvert(modVals,finalVals)
                                else newModVals <- .downConvert(modVals,finalVals)
                                #newModVals <- sort(newModVals)
                                #for (mvNum in seq(length(newModVals))) {
                                #    if ((newModVals[mvNum] <= lowVal[mvNum]) || (newModVals[mvNum] >= highVal[mvNum])) {
                                #        newModVals[mvNum] <- finalVals[mvNum]
                                #    }
                                #}
                                resListUpdate[[modName]] <- sort(newModVals)
                            }
                        }
                    }
                    #finally deal with expermential units with NA thresholds
                    #these are experimental units that did not produce any thresholds in the
                    #annotation forest
                    naNames <- intersect(names(which(is.na(resListUpdate))),allExpUnitsInIH)
                    if (length(naNames)) {
                        standardizedExpUnits <- setdiff(allExpUnitsInIH,naNames)
                        #recompute the gate matrix to account for expUnits that are now standard
                        gateMatrix <- matrix(nrow=0,ncol=gateNumber)
                        for (sExpUnit in standardizedExpUnits) {
                            gateData <- sort(resListUpdate[[sExpUnit]])
                            gateMatrix <- rbind(gateMatrix,gateData)
                            rownames(gateMatrix)[nrow(gateMatrix)] <- sExpUnit
                        }
                        finalVals <- apply(gateMatrix,2,stats::median)
                        for (changeName in naNames) {
                            resListUpdate[[changeName]] <- sort(finalVals)
                        }
                    }
                }
            }
            if (length(which(is.na(resListUpdate)))) {
                #for experimental units that still do not have
                #the standard gateNumber of annotation boundaries for the marker,
                #we now standardize across the entire experiment
                gateMatrix <- matrix(nrow=0,ncol=gateNumber)
                for (level in names(which(!is.na(resListUpdate)))) {
                    gateData <- sort(resListUpdate[[level]])
                    gateMatrix <- rbind(gateMatrix,gateData)
                    rownames(gateMatrix)[nrow(gateMatrix)] <- level
                }
                finalVals <- apply(gateMatrix,2,stats::median)
                mgMAD <- apply(gateMatrix,2,stats::mad)
                if (any(mgMAD == 0)) {
                    mgMAD <- Inf
                }
                lowVal <- finalVals - (stats::qt(0.975,df=1) * mgMAD)
                highVal <- finalVals + (stats::qt(0.975,df=1) * mgMAD)
                #map those still NA to the experiment wide medians.
                naNames <- names(which(is.na(resListUpdate)))
                for (changeName in naNames) {
                    modVals <- gateList[[changeName]][[channel]]
                    if (any(!is.na(modVals))) {
                        #there is empirical data for an experimental unit, so attempt to use it.
                        #this can arise if the "Preference" setting of supervision leads to using
                        #a gateNumber that is absent from a level of the imputation hierarchy.
                        if (gateNumber < length(modVals)) newModVals <- .upConvert(modVals,finalVals)
                        else newModVals <- .downConvert(modVals,finalVals)
                        newModVals <- sort(newModVals)
                        for (mvNum in seq(length(newModVals))) {
                            if ((newModVals[mvNum] <= lowVal[mvNum]) || (newModVals[mvNum] >= highVal[mvNum])) {
                                newModVals[mvNum] <- finalVals[mvNum]
                            }
                        }
                        resListUpdate[[changeName]] <- sort(newModVals)
                    }
                    else {
                        #impute the values outright
                        resListUpdate[[changeName]] <- sort(finalVals)
                    }
                }
            }
            resList[[channel]] <- resListUpdate
        }
        if (length(forceList) > 0) {
            #the user has indicated a channel must be included in the anlaysis and gated at a value.
            #add it in now, overwriting any automatic reconciliation.
            forcedNames <- names(forceList)
            listTemplate <- resList[[1]]
            designSettings <- names(listTemplate)
            selectedChannelsOut <- selectedChannels
            for (forcedMarkerName in forcedNames) {
                forcedGates <- forceList[[forcedMarkerName]] #the forced gate values
                newTemplate <- listTemplate #copy the template for updating
                for (setting in designSettings) {
                    newTemplate[[setting]] <- forcedGates
                }
                if (forcedMarkerName %in% selectedChannelsOut) {
                    print(paste0("Overwriting empirical gates for user settings on marker ",
                                 forcedMarkerName))
                    resList[[forcedMarkerName]] <- newTemplate
                }
                else {
                    resList <- append(resList,list(newTemplate))
                    names(resList)[length(resList)] <- forcedMarkerName
                    selectedChannelsOut <- append(selectedChannelsOut,
                                                  forcedMarkerName)
                }
            }
        }
        else {
            selectedChannelsOut <- selectedChannels
        }
        saveRDS(resList,
                file.path(normalizePath(projectPath),
                          "faustData",
                          "gateData",
                          paste0(parentNode,"_resListPrep.rds")))
        saveRDS(selectedChannelsOut,
                file.path(normalizePath(projectPath),
                          "faustData",
                          "gateData",
                          paste0(parentNode,"_selectedChannels.rds")))
    }
    return()
}

.getConsistentNumberOfGatesForMarkers <- function(numGateMatrix,preferenceList,projectPath,debugFlag)
{
    numSel <- c()
    possibilityList <- list()
    if (length(preferenceList) > 0) preferredMarkers <- names(preferenceList)
    for (columnName in colnames(numGateMatrix)) {
        nzLookup <- which(numGateMatrix[,columnName] > 0)
        if (length(nzLookup)==0) {
            print(numGateMatrix)
            print(columnName)
            stop("Invalid reconciliation")
        }
        columnSub <- numGateMatrix[nzLookup,columnName]
        columnCounts <- table(columnSub)
        possibilityList <- append(possibilityList,list(columnCounts))
        names(possibilityList)[length(possibilityList)] <- columnName
        columnMax <- max(columnCounts)
        if ((length(preferenceList) > 0) && (columnName %in% preferredMarkers)) {
            #user has requested a particular number of gates for this marker.
            #if there is empirical evidence for this choice, accommodate it.
            #otherwise, defer to max choice.
            columnPreference <- preferenceList[[columnName]]
            preferenceLookup <- which(as.numeric(names(columnCounts)) == as.numeric(columnPreference))
            if (length(preferenceLookup) > 0) {
                #we have observed the preference in the data.
                selUpdate <- as.numeric(columnPreference)
                if (debugFlag) print(paste0("Selecting ",selUpdate," gates for marker ",columnName))
            }
            else {
                #the request is not supported by the data. Use standard rule.
                selUpdate <- sort(as.numeric(names(which(columnCounts==columnMax))))[1]
                if (debugFlag) {
                    print(paste0("No empirical evidence for user gating preference in marker ",columnName))
                    print(paste0("Selecting ",selUpdate," gates for marker ",columnName))
                }
            }
        }
        else {
            selUpdate <- sort(as.numeric(names(which(columnCounts==columnMax))))[1] #if non-zero ties, pick min
        }
        numSel <- append(numSel,selUpdate)
        names(numSel)[length(numSel)] <- columnName
    }
    #save the possible number of gates for all selected markers in the experiment.
    #users can examine this data to override max selection
    saveRDS(possibilityList,
            file.path(normalizePath(projectPath),
                      "faustData",
                      "metaData",
                      "possibilityList.rds"))
    return(numSel)
}


.makeGateList <- function(uniqueExpUnits,projectPath,parentNode,selectedChannels)
{
    gateList <- list()
    for (expUnit in uniqueExpUnits) {
        if (file.exists(file.path(normalizePath(projectPath),
                                  "faustData",
                                  "expUnitData",
                                  expUnit,
                                  paste0(parentNode,"_pAnnF.rds"))))
        {
            afIn <- readRDS(file.path(normalizePath(projectPath),
                                      "faustData",
                                      "expUnitData",
                                      expUnit,
                                      paste0(parentNode,"_pAnnF.rds")))
            gates <- lapply(selectedChannels,
                            function(x){eval(parse(text=paste0("afIn$`",x,"`$gates")))})
            if (length(gates) == length(selectedChannels)) {
                names(gates) <- selectedChannels
                gateList <- append(gateList,list(gates))
                names(gateList)[length(gateList)] <- expUnit
            }
        }
        else {
            print(paste0(expUnit,": no parsed forest detected."))
        }
    }
    return(gateList)
}

.upConvert <- function(fromGates,toGates) {
    #toGates is assumed sorted in increasing order.
    dropInds <- c()
    outGates <- c()
    for (gate in fromGates) {
        #for each fromGate, find closest to gate and record index
        diffVec <- abs(toGates-gate)
        minDiff <- which(diffVec==min(diffVec))
        dropInds <- append(dropInds,minDiff[1]) #always keep higher gates.
        outGates <- append(outGates,gate)
    }
    #if multiple fromGates match a single toGate, use smallest fromGate
    uniqDrops <- sort(unique(dropInds))
    if (length(uniqDrops) != length(dropInds)) {
        finalOG <- c()
        for (indexVal in uniqDrops) {
            finalOG <- append(finalOG,outGates[which(dropInds == indexVal)][1])
        }
    }
    else {
        finalOG <- outGates
    }
    #copy over indices which are not explained
    #note length(uniqDrops)>0 since length(toGates) > length(fromGates)
    finalOG <- sort(unique(append(finalOG,toGates[-uniqDrops])))
    if (length(finalOG) != length(toGates)) {
        stop("Upconversion error.")
    }
    return(finalOG)
}

.downConvert <- function(fromGates,toGates) {
    #toGates is assumed sorted in increasing order.
    minInds <- c()
    minDiffVals <- c()
    for (gate in fromGates) {
        diffVec <- abs(toGates-gate)
        minVal <- min(diffVec)
        minDiff <- which(diffVec==minVal)
        minInds <- append(minInds,minDiff[1])#always match to lower gate in ties
        minDiffVals <- append(minDiffVals,minVal)
    }
    uniqMinInds <- sort(unique(minInds))
     #take the best of the fromGates that map to toGates
    finalOG <- c()
    for (indexVal in uniqMinInds) {
        subLookup <- which(minInds==indexVal)
        subDiffVals <- minDiffVals[subLookup]
        subGates <- fromGates[subLookup]
        allKeepInd <- which(subDiffVals == min(subDiffVals))
        keepInd <- allKeepInd[length(allKeepInd)]
        finalOG <- append(finalOG,subGates[keepInd])
    }
    #add on the toGates unaccounted for
    if (length(uniqMinInds) != length(toGates)) {
        remInds <- setdiff(seq(length(toGates)),uniqMinInds)
        finalOG <- append(finalOG,toGates[remInds])
    }
    finalOG <- sort(unique(finalOG))
    if (length(finalOG) != length(toGates)) {
        stop("Downconversion error.")
    }
    return(finalOG)
}
