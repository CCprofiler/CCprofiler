context("Consecutive id filter")


test_that("Fringe cases", {
  testthat::expect_error(filterConsecutiveIdStretches(c()),"Object is not of class traces.")

})
test_that("Output Format", {
  testPeptides <- examplePeptideTraces$traces$id[1:3]
  testTraces <- subset(examplePeptideTraces, trace_subset_ids = testPeptides)
  
  tracesFilteredConSec <- filterConsecutiveIdStretches(testTraces,
                                                       min_stretch_length = 3,
                                                       remove_empty = T)
  
  tracesFilteredConSecNoRM <- filterConsecutiveIdStretches(examplePeptideTraces,
                                                           min_stretch_length = 3,
                                                           remove_empty = F)
  testthat::expect_identical(tracesFilteredConSec$trace_annotation,
                             subset(examplePeptideTraces, testPeptides[c(1,3)])$trace_annotation)
  testthat::expect_identical(tracesFilteredConSecNoRM$trace_annotation,
                             examplePeptideTraces$trace_annotation)
  
})
