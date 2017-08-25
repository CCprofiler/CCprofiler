

test_that("imputeMissingValues",{
  tracesMatrix <- getIntensityMatrix(traces)
  tracesImputed_5p <- imputeMissingValues(tracesMatrix,"5%")
  tracesImputed_val <- imputeMissingValues(tracesMatrix,1587)
  testthat::expect_equal(tracesImputed_5p,tracesImputed_val)
  testthat::expect_equal(sum(tracesImputed_5p==0),0)
})
