Posdef <- function (n, ev = runif(n, 1, 2))
{
    #function written by Ravi Varadhan, from r-help mailing list
    #Thu Feb 7, 20:02:30 CET 2008
    #Generates a positive a positive definine matrix of dimension n
    #ev bounds the covariance by 2
    Z <- matrix(ncol=n, rnorm(n^2))
    decomp <- qr(Z)
    Q <- qr.Q(decomp)
    R <- qr.R(decomp)
    d <- diag(R)
    ph <- d / abs(d)
    O <- Q %*% diag(ph)
    Z <- t(O) %*% diag(ev) %*% O
    return(Z)
}

labelMeanVector <- function(meanVector) {
    baseStr <- paste0("V",seq(length(meanVector)))
    tokenStr <- rep("+",length(meanVector))
    tokenStr[which(meanVector==0)] <- "-"
    return(paste0(paste0(baseStr,tokenStr),collapse=""))
}

genClusterCentersWithLabels <- function(possibleCenterMat,
                                        nPop=5,
                                        seedVal=0)
{
    #given the possibleCenterMat, 
    #sample a location value from each column to determine 
    #the mean vector of a cluster.
    #annotations are then determined by the sampled meanvector
    #Returns the fixedMeanMatrix, which locates each cluster's center
    #Returns the fixedLabelVector, which is the true population of the cluster.
    pcMat <- t(possibleCenterMat)
    allPops <- expand.grid(split(pcMat,rep(seq(nrow(pcMat)),ncol(pcMat))))
    if (nPop > nrow(allPops)) {
        print("Too many clusters relative to specMat/possibleCenterMat.")
        stop("Reduce number of clusters or increase possible mean vectors.")
    }
    clusterIndex <- sample(seq(nrow(allPops)),nPop)
    outMat <- allPops[clusterIndex,,drop=FALSE]
    outLab <- as.character(apply(outMat,1,labelMeanVector))
    outList <- list(fixedMeanMatrix=outMat,fixedLabelVector=outLab)
    return(outList)
}

simSample <- function(sampleDim,
                      sampleSize,
                      batchEffect,
                      transformationType,
                      sRegime=0,
                      mixtureType="gaussianOnly",
                      fixedMeanMatrix=NA,
                      fixedLabelVector=NA,
                      noiseDim=0,
                      probVecSample,
                      targetRank)
{
    numClusters <- nrow(fixedMeanMatrix)
    probVec <- probVecSample
    if (sRegime) {
        #in the spiked in regime, we increase the prevalance of the spiked population.
        currentSpikeMass <- probVec[length(probVec)]
        targetSpikeMass <- sort(probVec,decreasing=TRUE)[targetRank]
        #spread existing mass equally across all populations.
        massIncrement <- currentSpikeMass/(length(clusterProbVec)-1)
        probVec <- probVec + massIncrement
        #zero out the spike
        probVec[length(probVec)] <- 0
        #proportionally decrement mass by the target
        probVec <- probVec - (targetSpikeMass * probVec)
        probVec[length(probVec)] <- targetSpikeMass
        if (abs(sum(probVec) - 1) > 1e-9) {
            print(probVec)
            print(sum(probVec))
            print(abs(sum(probVec) - 1))
            stop("Error in probability reapportionment")
        }
        if (min(probVec) == probVec[length(probVec)]) {
            stop("Error in spiking in population")
        }
    }
    sampleSizeVec <- as.vector(t(rmultinom(1,sampleSize,probVec)))
    outData <- matrix(nrow=0,ncol=(sampleDim+noiseDim))
    outLabel <- c()
    #for (sampleNum in seq(length(sampleSizeVec))) {
    for (sampleNum in seq(numClusters)) {
        #if we are simulating a sample without the spiked in pop, skip it.
        if (sampleSizeVec[sampleNum] == 0) {
            #print("skip")
            next
        }
        currentMu <- as.numeric(fixedMeanMatrix[sampleNum,,drop=TRUE])
        currentMu <- currentMu + batchEffect #a constant batch effect translates the sample mean.
        subjectShift <- round(rnorm(length(currentMu),mean=0,sd=(1/sqrt(2))))
        currentMu <- currentMu + subjectShift #model sample specifc translations of the mean vector
        currentLabel <- fixedLabelVector[sampleNum]
        outLabel <- append(outLabel,rep(currentLabel,sampleSizeVec[sampleNum]))
        currentSigma<- Posdef(sampleDim)
        if (mixtureType == "tPlusGauss") {
            if ((sampleNum %% 2) == 0) {
                currentSample <- mvrnorm(sampleSizeVec[sampleNum],mu=currentMu,Sigma=currentSigma)
            }
            else {
                currentSample <- rmvt(sampleSizeVec[sampleNum],delta=currentMu,sigma=currentSigma,df=5)
            }
        }
        else  {
            currentSample <- mvrnorm(sampleSizeVec[sampleNum],mu=currentMu,Sigma=currentSigma)
        }
        if (is.vector(currentSample)) {
            currentSample <- t(as.matrix(currentSample))
        }
        if (noiseDim > 0) {
            currentSigma<- Posdef(noiseDim)
            noiseSample <- mvrnorm(nrow(currentSample),mu=rep(5,noiseDim),Sigma=currentSigma)
            if (is.vector(noiseSample)) {
                noiseSample <- t(as.matrix(noiseSample))
            }
            currentSample <- cbind(currentSample,noiseSample)
        }
        outData <- rbind(outData,t(apply(currentSample,1,transformationType)))
    }
    colnames(outData) <- paste0("V",seq(ncol(outData)))
    outList <- list(sampleMatrix=outData,sampleLabels=outLabel)
    return(outList)
}


getStrOrder <- function(faustLabel) {
    orderVec <- c()
    while (nchar(faustLabel) > 0) {
        currentToken <- substr(faustLabel,1,1)
        if (currentToken=="V") {
            orderStr <- c()
        }
        else if (currentToken %in% c("-","M","+")) {
            orderVec <- append(orderVec,as.numeric(orderStr))
        }
        else {
            orderStr <- paste0(append(orderStr,currentToken),collapse="")
        }
        if (nchar(faustLabel) == 1){
            faustLabel <- ""
        }
        else {
            faustLabel <- substr(faustLabel,2,nchar(faustLabel))
        }
    }
    return(orderVec)
}

truth2faust <- function(trueLabel,faustOrder) {
    #reorder a true label by the depth-score faustOrder.
    tokenStr <- labelVec <- c()
    vctr <- 0
    while (nchar(trueLabel) > 0) {
        currentToken <- substr(trueLabel,1,1)
        if (currentToken=="V") {
            vctr <- vctr + 1
            if (vctr > 1) {
                labelVec <- append(labelVec,tokenStr)
                tokenStr <- c()
            }
        }
        tokenStr <- paste0(append(tokenStr,currentToken),collapse="")
        if (nchar(trueLabel) == 1){
            trueLabel <- ""
        }
        else {
            trueLabel <- substr(trueLabel,2,nchar(trueLabel))
        }
    }
    labelVec <- append(labelVec,tokenStr)
    return(paste0(labelVec[faustOrder],collapse=""))
}

safeMod <- function(dataSet) {
    out <- tryCatch(
    {
        m <- glmer(cbind(childCount,(parentCount-childCount)) ~ resType + (1|subjectID),
                   data=dataSet,family="binomial",
                   control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1e5),
                                        check.conv.singular = .makeCC(action = "warning",  tol = 1e-4)))
        rv <- c(coefficients(summary(m))[2,4])
        return(rv)
    },
    error=function(cond){
        message("Error!")
        message(cond)
        return(NA)
    },
    warning=function(cond){
        message("Warning!")
        message(cond)
        return(NA)
    },
    finally={
        message("done!")
    }
    )
    return(out)
}

simulateExperiment <- function(meanVectorBoundsMatrix,
                               numSamples=100,
                               randomSeed=0,
                               transformationList=list(function(x){return(x)},
                                                       function(x){return(x)},
                                                       function(x){return(x)}),
                               batchEffectShift=0.25,
                               noiseDimension=5,
                               probVecIn,
                               minSampleSize=5000,
                               tncp=10000,
                               targetRank)
{
    sampleSpecs <- genClusterCentersWithLabels(
        possibleCenterMat=meanVectorBoundsMatrix,
        nPop=length(probVecIn),
        seedVal=randomSeed
    )
    currentTransformation <- transformationList[[1]]
    #the regime determines if we spike in a population
    #regime==1 -> spike it in.
    regime <- rep(0,numSamples)
    regime[sample(seq(numSamples),(numSamples/2))] <- 1
    currentSample <- 1
    nextBatchEffect <- rep((-1*batchEffectShift),ncol(sampleSpecs$fixedMeanMatrix))
    labelsList <- regimeList <- flowList <- truthList <- list()
    while (currentSample <= numSamples) {
        #bound an experimental sample at 5000 observations
        nextSampleSize <- max(minSampleSize,round(rt(1,df=3,ncp=tncp)))
        nextRegime <- regime[currentSample]
        if (((currentSample - 1) %% 10) == 0) {
            nextBatchEffect <- nextBatchEffect + batchEffectShift
            #print(paste0("Batch: ",nextBatchEffect[1]))
            if ((currentSample - 1) > floor(numSamples/3)) {
                currentTransformation <- transformationList[[2]]
            }
            if ((currentSample - 1) > floor((2*(numSamples/3)))) {
                currentTransformation <- transformationList[[3]]
            }
        }
        sampleData <- simSample(
            sampleDim=ncol(sampleSpecs$fixedMeanMatrix),
            sampleSize=nextSampleSize,
            batchEffect=nextBatchEffect,
            sRegime=nextRegime,
            transformationType=currentTransformation,
            fixedMeanMatrix=sampleSpecs$fixedMeanMatrix,
            fixedLabelVector=sampleSpecs$fixedLabelVector,
            noiseDim=noiseDimension,
            probVecSample=probVecIn,
            targetRank=targetRank
        )
        if (currentSample < 10) {
            outName <- paste0("sample00",currentSample)
        }
        else if (currentSample < 100) {
            outName <- paste0("sample0",currentSample)
        }
        else {
            outName <- paste0("sample",currentSample)
        }
        ff <- flowFrame(sampleData$sampleMatrix)
        flowList <- append(flowList,ff)
        names(flowList)[length(flowList)] <- outName
        labelsList <- append(labelsList,list(sampleData$sampleLabels))
        names(labelsList)[length(labelsList)] <- outName
        truthList <- append(truthList,list(table(sampleData$sampleLabels)))
        names(truthList)[length(truthList)] <- outName
        regimeList <- append(regimeList,list(nextRegime))
        names(regimeList)[length(regimeList)] <- outName
        currentSample <- currentSample + 1
    }
    truthMat <- getSimTrueCountMatrix(truthList)
    
    outputList <- list("flowFrameList"=flowList,
                       "truthValueList"=truthList,
                       "truthMat"=truthMat,
                       "spikedPop"=sampleSpecs$fixedLabelVector[length(probVecIn)],
                       "spikedInPopList"=regimeList,
                       "labelsList"=labelsList)
    return(outputList)
}

getSimTrueCountMatrix <- function(truthList)
{
    uniqueNames <- c()
    for (i in seq(length(truthList))) {
        uniqueNames <- append(uniqueNames,names(truthList[[i]]))
    }
    uniqueNames <- sort(unique(uniqueNames))
    truthMat <- matrix(0,nrow=length(truthList),ncol=length(uniqueNames))
    colnames(truthMat) <- uniqueNames
    rownames(truthMat) <- names(truthList)
    for (i in seq(length(truthList))) {
        cv <- truthList[[i]]
        for (cName in colnames(truthMat)) {
            lookup <- which(names(cv) == cName)
            if (length(lookup)) {
                truthMat[i,cName] <- cv[lookup]
            }
        }
    }
    return(truthMat)
}



#modified version of function "external_validation" in ClusterR package 
#https://cran.r-project.org/web/packages/ClusterR/index.html
#adjusted to provide methodName parameter
#this version is not affiliated with the ClusterR package or its auther/maintainer.
myExternalValidation = function(true_labels, clusters, methodName) {
    #modified from ClusterR package so it returns a column vector of clustering methods
    if (is.integer(true_labels)) true_labels = as.numeric(true_labels)
    if (is.integer(clusters)) clusters = as.numeric(clusters)
    if (!is.vector(true_labels) || !is.numeric(true_labels)) stop('true_labels should be a numeric vector')
    if (!is.vector(clusters) || !is.numeric(clusters)) stop('clusters should be a numeric vector')
    if (length(true_labels) != length(clusters)) stop('the length of the true_labels vector should equal the length of the clusters vector')

    tbl = table(clusters, true_labels)

    conv_df = as.data.frame.matrix(tbl)

                                        # Diagonal = rep(0, ncol(conv_df))
                                        #
                                        # for (i in 1:nrow(conv_df)) {
                                        #
                                        #   wh_idx = which.max(conv_df[i, ])
                                        #
                                        #   if (conv_df[i, wh_idx] > Diagonal[wh_idx]) {
                                        #
                                        #     Diagonal[wh_idx] = conv_df[i, wh_idx]
                                        #   }
                                        # }

    conv_df = as.data.frame.matrix(tbl)

    tp_plus_fp = sum(gmp::asNumeric(gmp::chooseZ(rowSums(conv_df), 2)))

    tp_plus_fn = sum(gmp::asNumeric(gmp::chooseZ(colSums(conv_df), 2)))

    tp = sum(gmp::asNumeric(gmp::chooseZ(as.vector(as.matrix(conv_df)), 2)))

    fp = tp_plus_fp - tp

    fn = tp_plus_fn - tp

    tn = gmp::asNumeric(gmp::chooseZ(sum(as.vector(as.matrix(conv_df))), 2)) - tp - fp - fn


    prod_comb = (tp_plus_fp * tp_plus_fn) / gmp::asNumeric(gmp::chooseZ(length(true_labels), 2))

    mean_comb = (tp_plus_fp + tp_plus_fn) / 2.0



    tmp_pur = apply(conv_df, 1, max)

    res_purity = sum(tmp_pur)/length(true_labels)


    tmp_entropy = sum(apply(conv_df, 2, function(x) entropyFormula(x)))

    res_entropy = -(1/(sum(tbl) * log2(length(unique(true_labels))))) * tmp_entropy


    mutual_information = 0.0

    joint_entropy = 0.0

    for (i in 1:nrow(conv_df)) {

        for (j in 1:ncol(conv_df)) {

            if (conv_df[i,j] > 0.0) {

                joint_entropy = joint_entropy + (-((conv_df[i,j] / sum(tbl)) * log2(conv_df[i,j] / sum(tbl))))

                mutual_information = mutual_information + ((conv_df[i,j] / sum(tbl)) * log2((sum(tbl) * conv_df[i,j]) / (sum(conv_df[i,]) * sum(conv_df[,j]))))
            }
        }
    }

    entr_cluster = sum(apply(conv_df, 1, function(x) -(sum(x) / sum(tbl)) * log2(sum(x) / sum(tbl))))

    entr_class = sum(apply(conv_df, 2, function(x) -(sum(x) / sum(tbl)) * log2(sum(x) / sum(tbl))))

    NMI = (mutual_information / ((entr_cluster + entr_class) / 2.0))

    VAR_INFO = (entr_cluster + entr_class) - 2.0 * mutual_information

    NVI = 1.0 - (mutual_information / joint_entropy)


    prec = tp / (tp + fp)
    rec = tp / (tp + fn)

    ovec <- c(round(res_purity, 3),
              round(res_entropy,3),
              round(NMI,3),
              round(VAR_INFO,3),
              round(NVI,3),
              round(tn / (tn + fp), 3),
              round(tp / (tp + fn), 3),
              round(prec,3),
              round(rec,3),
              round(2.0 * ((prec * rec) / (prec + rec)), 3),
              round((tp + tn) / (tp + fp + fn + tn), 3),
              round((tp - prod_comb) / (mean_comb - prod_comb), 3),
              round(tp / (tp + fp + fn), 3),
              round(sqrt((tp / ((tp + fp))) * (tp / (tp + fn))), 3))#,
              #round(2.0 * (fp + fn), 4))
    omat <- matrix(ovec,ncol=1)
    rownames(omat) <- c("Purity","Entropy","NormMutInfo","VarInfo","NormVarInfo","Specificity","Sensitivity",
                        "Precision","Recall","F-measure","Rand Index", "Adj Rand Index", "Jaccard Index",
                        "Fowlkes-Mallows Index")#, "Mirkin Metric")
    
    colnames(omat) <- methodName
    return(omat)
    
}


entropyFormula = function(x_vec) {
    #copied from ClusterR package. Used in modified external_validation call.
    vec = rep(NA, length(x_vec))
    for (i in 1:length(x <- vec)) {
        if (x_vec[i] == 0.0) {
            vec[i] = 0.0}
        else {
            vec[i] = ((x_vec[i]) * log2(x_vec[i]/sum(x_vec)))
        }
    }
    return(vec)
}
