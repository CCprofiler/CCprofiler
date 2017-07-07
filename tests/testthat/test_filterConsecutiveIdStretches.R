context("Consecutive id filter")

testTraces <- subset(examplePeptideTraces, c("AAAAPGASPSPGGDAAWSEAGPGPR",
                                             "AAAAPGASPSPGGDAAWSEAGPGPRPLAR"))

tracesFilteredConSec <- filterConsecutiveIdStretches(testTraces,
                                                     min_stretch_length = 3,
                                                     remove_empty = T)

tracesFilteredConSecNoRM <- filterConsecutiveIdStretches(examplePeptideTraces,
                                                     min_stretch_length = 3,
                                                     remove_empty = F)

test_that("Fringe cases", {
  testthat::expect_error(filterConsecutiveIdStretches(c()),"Object is not of class traces.")

})
test_that("Output Format", {
  testthat::expect_identical(tracesFilteredConSec$trace_annotation,
                             subset(examplePeptideTraces, "AAAAPGASPSPGGDAAWSEAGPGPRPLAR")$trace_annotation)
  testthat::expect_identical(tracesFilteredConSecNoRM$trace_annotation,
                             examplePeptideTraces$trace_annotation)
  
})
