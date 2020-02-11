context("Test functionality related to standardizing annotation boundaries.")

test_that("Up conversion works for multiple different standards.",{
    expect_equal(faust:::.upConvert(c(5),c(2.5,8)),c(5,8))
    expect_equal(faust:::.upConvert(c(5),c(2.5,7.5,10)), c(5,7.5,10))
    expect_equal(faust:::.upConvert(c(5,6),c(2.5,10,15)), c(5,10,15))
    expect_equal(faust:::.upConvert(c(5,6),c(2.5,8,15)), c(5,6,15))
})

test_that("Down conversion down coversion works for multiple different standards.",{
    expect_equal(.downConvert(c(1,2,3,4),c(5,10)),c(4,10))
    expect_equal(.downConvert(c(1,2,3,4),c(0,5,10)),c(1,4,10))
    expect_equal(.downConvert(c(-1,1,2,3,4),c(0,5,10)),c(1,4,10))
})
