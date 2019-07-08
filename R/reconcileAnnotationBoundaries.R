.reconcileAnnotationBoundaries <- function(selectedChannels,
                                           parentNode,
                                           analysisMap,
                                           initSplitPval,
                                           projectPath,
                                           debugFlag,
                                           preferenceList,
                                           forceList,
                                           manualList) {
    if (!dir.exists(paste0(projectPath,"/faustData/gateData"))) {
        dir.create(paste0(projectPath,"/faustData/gateData"))
    }
    uniqueLevels <- unique(analysisMap[,"analysisLevel"])
    gateList <- .makeGateList(
        uniqueLevels=uniqueLevels,
        projectPath=projectPath,
        parentNode=parentNode,
        selectedChannels=selectedChannels
    )
    if (length(gateList)) {
        saveRDS(gateList,paste0(projectPath,"/faustData/gateData/",parentNode,"_rawGateList.rds"))
        numGateMatrix <- Reduce(rbind,lapply(gateList,
                                             function(x){unlist(lapply(x,
                                                                       function(y){ifelse(is.na(y[1]),0,length(y))}))}))
        if (!is.matrix(numGateMatrix)) numGateMatrix <- t(as.matrix(numGateMatrix))
        rownames(numGateMatrix) <- names(gateList)
        numSel <- .getConsistentNumberOfGatesForMarkers(
            numGateMatrix=numGateMatrix,
            preferenceList=preferenceList,
            projectPath=projectPath,
            debugFlag=debugFlag
        )
        resList <- rep(list(NA),ncol(numGateMatrix))
        names(resList) <- colnames(numGateMatrix)
        for (channel in names(numSel)) {
            resListUpdate <- rep(list(NA),nrow(numGateMatrix))
            names(resListUpdate) <- rownames(numGateMatrix)
            #get data for our selected number of gates
            gateNumber <- as.numeric(numSel[channel])
            numLookup <- which(numGateMatrix[,channel]==gateNumber)
            matchLevels <- rownames(numGateMatrix)[numLookup]
            gateMatrix <- matrix(nrow=0,ncol=gateNumber)
            for (level in matchLevels) {
                gateData <- sort(gateList[[level]][[channel]])
                gateMatrix <- rbind(gateMatrix,gateData)
                rownames(gateMatrix)[nrow(gateMatrix)] <- level
            }
            fragmentedGates <- FALSE
            possibleGates <- as.numeric(names(table(numGateMatrix[,channel])))
            possibleGates <- setdiff(possibleGates,c(0,gateNumber))
            if (initSplitPval > 0)
            {
                #user has indicated they wish to test gates for homogeneity
                stillChecking <- TRUE
                numSplits <- 1
                while (stillChecking) {
                    anySplit <- FALSE
                    for (colNum in seq(ncol(gateMatrix))) {
                        naLookup <- which(is.na(gateMatrix[,colNum]))
                        if (length(naLookup)) dipPval <- singleDip(sort(gateMatrix[-naLookup,colNum]))
                        else dipPval <- singleDip(sort(gateMatrix[,colNum]))
                        if ((dipPval <= (initSplitPval/numSplits)))
                        {
                            numSplits <- numSplits + 1 
                            if (length(naLookup)) splitGate <- tsGates(sort(gateMatrix[-naLookup,colNum]),2)[2]
                            else splitGate <- tsGates(sort(gateMatrix[,colNum]),2)[2]
                            updateMatrix <- matrix(NA,nrow=nrow(gateMatrix),ncol=(ncol(gateMatrix)+1))
                            rownames(updateMatrix) <- rownames(gateMatrix)
                            for (upColNum in setdiff(seq(ncol(updateMatrix)),c(colNum,(colNum+1))))
                            {
                                if (upColNum < colNum) updateMatrix[,upColNum] <- gateMatrix[,upColNum]
                                else updateMatrix[,upColNum] <- gateMatrix[,(upColNum-1)]
                            }
                            bigUpdate <- gateMatrix[which(gateMatrix[,colNum] >= splitGate),colNum]
                            smallUpdate <- gateMatrix[which(gateMatrix[,colNum] < splitGate),colNum]
                            updateMatrix[which(gateMatrix[,colNum] >= splitGate),(colNum+1)] <- bigUpdate
                            updateMatrix[which(gateMatrix[,colNum] < splitGate),colNum] <- smallUpdate
                            newNumLookup <- which(numGateMatrix[,channel]==ncol(updateMatrix))
                            if (length(newNumLookup)) {
                                matchLevels <- rownames(numGateMatrix)[newNumLookup]
                                newNumMatrix <- matrix(nrow=0,ncol=ncol(updateMatrix))
                                for (level in matchLevels) {
                                    gateData <- sort(gateList[[level]][[channel]])
                                    newNumMatrix <- rbind(newNumMatrix,gateData)
                                    rownames(newNumMatrix)[nrow(newNumMatrix)] <- level
                                }
                                updateMatrix <- rbind(updateMatrix,newNumMatrix)
                            }
                            possibleGates <- setdiff(possibleGates,ncol(updateMatrix))
                            gateMatrix <- updateMatrix
                            fragmentedGates <- TRUE
                            anySplit <- TRUE
                            print(channel)
                            break
                        }
                    }
                    if (!anySplit) {
                        stillChecking <- FALSE
                    }
                }
            }
            if (fragmentedGates) {
                #if we did split our selection in the test for homogeneity, overwrite NA with gate median
                medVals <- apply(updateMatrix,2,function(x){stats::median(x,na.rm=TRUE)})
                for (colNum in seq(ncol(updateMatrix))) {
                    naLookup <- which(is.na(updateMatrix[,colNum]))
                    if (length(naLookup)) updateMatrix[naLookup,colNum] <- medVals[colNum]
                }
                updateMatrix[which(is.na(updateMatrix[,2])),2] <- medVals[2]
                gateMatrix <- updateMatrix
            }
            #update the resList
            for (level in rownames(gateMatrix)) {
                resListUpdate[[level]] <- sort(gateMatrix[which(rownames(gateMatrix)==level),])
            }
            #check for outliers
            mgMed <- apply(gateMatrix,2,stats::median)
            if (fragmentedGates) mgMAD <- apply(gateMatrix,2,stats::sd)
            else mgMAD <- apply(gateMatrix,2,stats::mad)
            #if we set the gate using only one example (via supervision), the MAD is 0. 
            #set to (-Inf,Inf) because want to keep whats found by annotation forest. 
            if (mgMAD == 0) mgMAD <- Inf
            lowVal <- mgMed - (stats::qt(0.975,df=1) * mgMAD)
            highVal <- mgMed + (stats::qt(0.975,df=1) * mgMAD)
            chkMatrix <- Reduce(cbind,lapply(seq(ncol(gateMatrix)),
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
            #for levels of analysis with different numbers of gates than the selection,
            #map them to selected values.
            if (length(possibleGates)) {
                for (gateNum in possibleGates) {
                    modLookup <- which(numGateMatrix[,channel]==gateNum)
                    modSamples <- rownames(numGateMatrix)[modLookup]
                    for (modName in modSamples) {
                        modVals <- gateList[[modName]][[channel]]
                        if (gateNum < length(finalVals)) newModVals <- .upConvert(modVals,finalVals)
                        else newModVals <- .downConvert(modVals,finalVals)
                        newModVals <- sort(newModVals)
                        for (mvNum in seq(length(newModVals))) {
                            if ((newModVals[mvNum] <= lowVal[mvNum]) || (newModVals[mvNum] >= highVal[mvNum])) {
                                newModVals[mvNum] <- finalVals[mvNum]
                            }
                        }
                        resListUpdate[[modName]] <- sort(newModVals)
                    }
                }
            }
            #finally deal with NA gates
            naNames <- names(which(is.na(resListUpdate)))
            if (length(naNames)) {
                for (changeName in naNames) {
                    resListUpdate[[changeName]] <- sort(finalVals)
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
            selectedChannelsPrep <- selectedChannels
            for (forcedMarkerName in forcedNames) {
                forcedGates <- forceList[[forcedMarkerName]] #the forced gate values
                newTemplate <- listTemplate #copy the template for updating
                for (setting in designSettings) {
                    newTemplate[[setting]] <- forcedGates
                }
                if (forcedMarkerName %in% selectedChannelsPrep) {
                    print(paste0("Overwriting empirical gates for user settings on marker ",
                                 forcedMarkerName))
                    resList[[forcedMarkerName]] <- newTemplate
                }
                else {
                    resList <- append(resList,list(newTemplate))
                    names(resList)[length(resList)] <- forcedMarkerName
                    selectedChannelsPrep <- append(selectedChannelsPrep,
                                                  forcedMarkerName)
                }
            }
        }
        else {
            selectedChannelsPrep <- selectedChannels
        }
        if (length(manualList) > 0) {
            #the user has provided manual gates for some marker.
            #add these gates in now.
            #one of two conditions must be met.
            #either the marker must be selected automatically by FAUST.
            #or manual gates must be provided for ALL levels by the user.
            #if the case the marker is selected by FAUST, a subsets of levels may be manually gated.
            manualNames <- names(manualList)
            listTemplate <- resList[[1]]
            designSettings <- names(listTemplate)
            selectedChannelsOut <- selectedChannelsPrep
            for (manualMarkerName in manualNames) {
                manualGateList <- manualList[[manualMarkerName]] #the forced gate values
                manualSettings <- names(manualGateList)
                if (manualMarkerName %in% selectedChannelsOut) {
                    #the marker is selected and gated by FAUST. 
                    #get faust gates and update only those levels modified by the user.
                    newTemplate <- resList[[manualMarkerName]]
                    for (setting in manualSettings) {
                        newTemplate[[setting]] <- manualGateList[[setting]]
                    }
                    resList[[manualMarkerName]] <- newTemplate
                }
                else {
                    #the marker is not selected by FAUST
                    #therefore, the user must set gates for all design levels.
                    #check to make sure this occurs before updating.
                    newTemplate <- listTemplate
                    if (length(setdiff(designSettings,manualSettings))) {
                        print("Manual gates do not include all levels.")
                        print("Set gates for the following.")
                        print(setdiff(designSettings,manualSettings))
                        stop("Killing faust run.")
                    }
                    if (length(setdiff(manualSettings,designSettings))) {
                        print("Manual gates include levels not specified by the design.")
                        print("These levels are the following.")
                        print(setdiff(manualSettings,designSettings))
                        stop("Killing faust run.")
                    }
                    for (setting in manualSettings) {
                        newTemplate[[setting]] <- manualGateList[[setting]]
                    }
                    resList <- append(resList,list(newTemplate))
                    names(resList)[length(resList)] <- manualMarkerName
                    selectedChannelsOut <- append(selectedChannelsOut,
                                                  manualMarkerName)
                }
            }
        }
        else {
            selectedChannelsOut <- selectedChannelsPrep
        }
        saveRDS(resList,paste0(projectPath,"/faustData/gateData/",parentNode,"_resListPrep.rds"))
        saveRDS(selectedChannelsOut,paste0(projectPath,"/faustData/gateData/",parentNode,"_selectedChannels.rds"))
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
    saveRDS(possibilityList,paste0(projectPath,"/faustData/metaData/possibilityList.rds"))
    return(numSel)
}


.makeGateList <- function(uniqueLevels,projectPath,parentNode,selectedChannels)
{
    gateList <- list()
    for (aLevel in uniqueLevels) {
        if (file.exists(paste0(projectPath,"/faustData/levelData/",aLevel,"/",parentNode,"_pAnnF.rds"))) {
            afIn <- readRDS(paste0(projectPath,"/faustData/levelData/",aLevel,"/",parentNode,"_pAnnF.rds"))
            gates <- lapply(selectedChannels,
                            function(x){eval(parse(text=paste0("afIn$`",x,"`$gates")))})
            if (length(gates) == length(selectedChannels)) {
                names(gates) <- selectedChannels
                gateList <- append(gateList,list(gates))
                names(gateList)[length(gateList)] <- aLevel
            }
        }
        else {
            print(paste0(aLevel,": no parsed forest detected."))
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
#.upConvert(c(5),c(2.5,8)) == c(5,8)
#.upConvert(c(5),c(2.5,7.5,10)) == c(5,7.5,10)
#.upConvert(c(5,6),c(2.5,10,15)) == c(5,10,15)
#.upConvert(c(5,6),c(2.5,8,15)) == c(5,6,15)

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
#.downConvert(c(1,2,3,4),c(5,10)) == c(4,10)
#.downConvert(c(1,2,3,4),c(0,5,10)) == c(1,4,10)
#.downConvert(c(-1,1,2,3,4),c(0,5,10)) == c(1,4,10)
