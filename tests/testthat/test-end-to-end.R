context("end-to-end")

test_that("faust runs end to end with multiple samples.", {
  # Generate some fake data
  library(flowCore)
  library(flowWorkspace)
  means <- c(1,5,1,1,5,5,1)
  mats <- replicate(2,matrix(rnorm(n = 100000, mean = means, sd = 0.75), ncol = 2, byrow = TRUE))
  colnames(mats) <- c("A","B")
  f <- sapply(1:(dim(mats)[3]),function(i)flowFrame(mats[,,i]))
  
  gs <- GatingSet(flowSet(f))
  unlink(x = file.path(tempdir(),"faustData"),recursive = TRUE, force = TRUE)
  #test faust end to end
  faust.output <- capture.output(expect_null(faust(
    gatingSet = gs,
    experimentalUnit = "name",
    activeChannels = c("A", "B"),
    startingCellPop = "root",
    projectPath = tempdir(),
    depthScoreThreshold = 0.01,
    selectionQuantile = 1.0,
    threadNum = 1,
    seedValue = 1234,
    debugFlag = FALSE,
    annotationsApproved = TRUE
  )))
  expect_true(file_test("-d",file.path(tempdir(),"faustData")))
  expect_true(file_test("-f",file.path(tempdir(),"faustData","faustCountMatrix.rds")))
  countMatrix <- readRDS(file.path(tempdir(),"faustData","faustCountMatrix.rds"))
  # test the clustering meets some expectations
  expect_true(all(colSums(countMatrix) > c(20000,20000,20000,10000,-1)))
})

test_that("faust runs on a single sample.", {
  library(flowCore)
  library(flowWorkspace)
  means <- c(1,5,1,1,5,5,1)
  mats <- replicate(1,matrix(rnorm(n = 1000, mean = means, sd = 0.75), ncol = 2, byrow = TRUE))
  colnames(mats) <- c("A","B")
  f <- sapply(1:(dim(mats)[3]),function(i)flowFrame(mats[,,i]))
  
  gs <- GatingSet(flowSet(f))
  unlink(x = file.path(tempdir(),"faustData"),recursive = TRUE, force = TRUE)
  #test faust end to end
  faust.output <- capture.output(expect_null(faust(
    gatingSet = gs,
    experimentalUnit = "name",
    activeChannels = c("A", "B"),
    startingCellPop = "root",
    projectPath = tempdir(),
    depthScoreThreshold = 0.01,
    selectionQuantile = 1.0,
    threadNum = 1,
    seedValue = 1234,
    debugFlag = FALSE,
    annotationsApproved = TRUE
  )))
  # read the counts (TODO should test that it exists)
  expect_true(file_test("-d",file.path(tempdir(),"faustData")))
  expect_true(file_test("-f",file.path(tempdir(),"faustData","faustCountMatrix.rds")))
  countMatrix <- readRDS(file.path(tempdir(),"faustData","faustCountMatrix.rds"))
  # test the clustering meets some expectations
  expect_true(all(colSums(countMatrix) > c(100,100,100,50,-1)))
})
