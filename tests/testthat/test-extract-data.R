context("Test extracting data from gating set.")

test_that("data extracted by faust matches counts from gating set.",{
    suppressWarnings(suppressMessages(library(flowCore)))
    suppressWarnings(suppressMessages(library(flowWorkspace)))
    suppressWarnings(suppressMessages(library(faust)))
    means <- c(1,5,1,1,5,5,1)
    sizes <- c(3*rnbinom(n=5, size=1000, prob=0.05))
    mats <- lapply(sizes,
                   function(x){flowCore::flowFrame(matrix(rnorm(n=x,mean=means,sd=0.75),
                                                          ncol=3,
                                                          byrow=TRUE,
                                                          dimnames=list(c(),c("A","B","C"))))})
    gs <- GatingSet(as(mats,"flowSet"))
    unlink(x = file.path(tempdir(),"faustData"),recursive = TRUE, force = TRUE)
    testingProjPath <- normalizePath(tempdir())
    if (dir.exists(testingProjPath)) {
        dir.create(file.path(testingProjPath,"faustData"))
        dir.create(file.path(testingProjPath,"faustData","metaData"))
        dir.create(file.path(testingProjPath,"faustData","sampleData"))
        nullReturn <- faust:::.extractDataFromGS(
                                  gs=gs,
                                  activeChannels=c("A","B","C"),
                                  startingCellPop="root",
                                  projectPath=testingProjPath,
                                  debugFlag=FALSE
                              )
        expectedCounts <- rep(NA,length(gs))
        names(expectedCounts) <- sampleNames(gs)
        for (sn in names(expectedCounts)) {
            expectedCounts[sn] <- nrow(flowWorkspace::gh_pop_get_data(gs[sn],"root"))
        }
        observedCounts <- rep(NA,length(expectedCounts))
        names(observedCounts) <- names(expectedCounts)
        for (sn in names(observedCounts)) {
            observedCounts[sn] <- nrow(readRDS(file.path(testingProjPath,
                                                         "faustData",
                                                         "sampleData",
                                                         sn,
                                                         "exprsMat.rds")))
        }
    }
    unlink(x = file.path(testingProjPath),recursive = TRUE, force = TRUE)
    expect_identical(observedCounts,expectedCounts)
})
