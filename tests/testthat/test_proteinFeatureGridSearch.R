context("proteinFeatureGridSearch")

peptideTraces <- subset(examplePeptideTracesFiltered,
                        trace_subset_ids = unique(examplePeptideTraces$trace_annotation$protein_id)[3],
                        trace_subset_type = "protein_id")

## Perform a small grid search for 1 parameter combinations
# Depending on the computational resources this can take several minutes
gridList <- performProteinGridSearch(traces = peptideTraces,
                                     corrs = c(0.5, 0.95),
                                     windows = 12,
                                     smoothing = 9,
                                     rt_heights = 5,
                                     n_cores = 2)

singleSearch <- findProteinFeatures(peptideTraces, perturb_cutoff = "5%")

test_that("performProteinGridSearch",{
  
  testthat::expect_error(performProteinGridSearch(c()), "Object is not of class traces.")
  
  testthat::expect_equal(length(gridList), 2)
  testthat::expect_identical(names(gridList[[1]]), names(gridList[[2]]))
  # @TODO these should eventually work
  # # These colnames are not conserved due to parallelization
  # colnames <- names(singleSearch)[!names(singleSearch) %in% c("area", "peak_corr")]
  # testthat::expect_identical(gridList[[2]][,colnames, with = F], singleSearch[,colnames, with = F])
  # testthat::expect_identical(gridList[[2]][,colnames, with = F], features[,colnames, with = F])
  # testthat::expect_identical(gridList[[2]][,colnames, with = F],
  #                            subset(exampleProteinFeatures, grepl(unique(peptideTraces$trace_annotation$protein_id),
  #                                                                exampleProteinFeatures$protein_id))[,colnames, with = F])

})

test_that("runGridProteinFeatureFinding",{
  params <- c(corr = 0.95, window = 12, smoothing = 9, rt_height = 5)
  features <- .runGridProteinFeatureFinding(params, pepTraces = peptideTraces)
  testthat::expect_identical(features[, names(features) %in% names(singleSearch), with = F],
                             singleSearch)
  testthat::expect_equal(ncol(features) - length(params), ncol(singleSearch))
  testthat::expect_true(all(features$corr == params["corr"]))
  testthat::expect_true(all(features$window == params["window"]))
  testthat::expect_true(all(features$smoothing_length == params["smoothing"]))
  testthat::expect_true(all(features$rt_height == params["rt_height"]))
  
})

