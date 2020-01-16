.getFaustCountMatrix <- function(analysisMap,
                                 selectedChannels,
                                 debugFlag=FALSE,
                                 projectPath=".")
{
    faustClusterNames <- readRDS(file.path(normalizePath(projectPath),
                                           "faustData",
                                           "metaData",
                                           "scampClusterNames.rds"))
    faustClusterNames <- append(faustClusterNames,"0_0_0_0_0")
    activeSamples <- analysisMap[,"sampleName"]
    faustCountMatrix <- matrix(0,nrow=length(activeSamples),ncol=length(faustClusterNames))
    rownames(faustCountMatrix) <- activeSamples
    colnames(faustCountMatrix) <- faustClusterNames
    for (sampleName in activeSamples) {
        if (debugFlag) print(paste0("Getting counts from sample ",sampleName))
        sNum <- which(rownames(faustCountMatrix)==sampleName)
        sAnn <- utils::read.table(file=file.path(normalizePath(projectPath),
                                                 "faustData",
                                                 "sampleData",
                                                 sampleName,
                                                 "faustAnnotation.csv"),
                                  header=F,
                                  sep="`",
                                  stringsAsFactors=FALSE)[,1]
        for (colName in colnames(faustCountMatrix)) {
            faustCountMatrix[sNum,colName] <- length(which(sAnn==colName))
        }
    }
    #pretty print the column names
    #save the original encoding for ploting gating strategies.
    newColNames <- faustColNames <- colnames(faustCountMatrix)
    newColNames <- as.character(sapply(newColNames,function(x){gsub("~1~[[:digit:]]~","-",x)}))
    newColNames <- as.character(sapply(newColNames,function(x){gsub("~2~2~","+",x)}))
    newColNames <- as.character(sapply(newColNames,function(x){gsub("~3~3~","Bright",x)}))
    newColNames <- as.character(sapply(newColNames,function(x){gsub("~4~4~","++",x)}))
    newColNames <- as.character(sapply(newColNames,function(x){gsub("~2~3~","Dim",x)}))
    newColNames <- as.character(sapply(newColNames,function(x){gsub("~2~4~","Med-",x)}))
    newColNames <- as.character(sapply(newColNames,function(x){gsub("~3~4~","Med+",x)}))
    colNameMap <- data.frame(faustColNames=faustColNames,
                             newColNames=newColNames,
                             stringsAsFactors=FALSE)
    colnames(faustCountMatrix) <- newColNames
    saveRDS(colNameMap,
            file.path(normalizePath(projectPath),
                      "faustData",
                      "metaData",
                      "colNameMap.rds"))
    saveRDS(faustCountMatrix,
            file.path(normalizePath(projectPath),
                      "faustData",
                      "faustCountMatrix.rds"))
    return()
}
