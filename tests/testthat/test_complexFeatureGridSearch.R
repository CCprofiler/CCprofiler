context("proteinFeatureGridSearch")

complexHypotheses <- subset(exampleComplexHypotheses, 
                            complex_id == exampleComplexHypotheses$complex_id[5])

## Perform a small grid search for 2 parameter combinations
gridList <- performComplexGridSearch(traces = exampleProteinTraces,
                                     complex_hypothesis = complexHypotheses,
                                     corrs = c(0.5, 0.95),
                                     windows = 15,
                                     smoothing = 11,
                                     rt_heights = 5,
                                     n_cores = 2)

singleSearch <- findComplexFeatures(exampleProteinTraces, complexHypotheses,
                                    corr_cutoff = 0.95, perturb_cutoff="5%", n_cores = 1,
                                    window_size = 15, smoothing_length = 11,
                                    rt_height = 5)


test_that("performProteinGridSearch",{

  testthat::expect_error(performComplexGridSearch(c()), 'Object is not of class traces.')
  testthat::expect_error(performComplexGridSearch(c(), c()), "Object is not of class traces.")

  testthat::expect_equal(length(gridList), 2)
  testthat::expect_identical(names(gridList[[1]]), names(gridList[[2]]))
  testthat::expect_true(all(gridList[[2]]$complex_id %in% gridList[[1]]$complex_id))
  
  # These colnames are not conserved due to parallelization
  # colnames <- names(singleSearch)[!names(singleSearch) %in% c("area", "peak_corr")] 
  # testthat::expect_identical(gridList[[2]][,colnames, with = F], singleSearch[,colnames, with = F])
  # testthat::expect_identical(gridList[[2]][,colnames, with = F],
  #                            subset(excomtr, complex_id == complexHypotheses$complex_id)[,colnames, with = F])
  
})

test_that("runGridComplexFeatureFinding",{
  params <- c(corr = 0.95, window = 15, smoothing = 11, rt_height = 5)
  features <- .runGridComplexFeatureFinding(params,
                                            protTraces = exampleProteinTraces,
                                            complex_hypothesis = complexHypotheses)
  testthat::expect_identical(features[, names(features) %in% names(singleSearch), with = F],
                             singleSearch)
  testthat::expect_equal(ncol(features) - length(params), ncol(singleSearch))
  testthat::expect_true(all(features$corr == params["corr"]))
  testthat::expect_true(all(features$window == params["window"]))
  testthat::expect_true(all(features$smoothing_length == params["smoothing"]))
  testthat::expect_true(all(features$rt_height == params["rt_height"]))
  
})
