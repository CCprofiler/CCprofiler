context("calibrateMW")

calibrationFunctions <- calibrateMW(exampleCalibrationTable,
                                       plot = FALSE,
                                       PDF = FALSE)

test_that("Fringe cases", {
  testthat::expect_error(calibrateMW(c()),"calibration_table input is neither file name or data.table")
  testthat::expect_error(calibrateMW(data.table()),"calibration_table seems to be empty")
  
})
test_that("Output Format", {
  testthat::expect_equal(calibrationFunctions$MWtoFraction(0), Inf)
  testthat::expect_lt(calibrationFunctions$MWtoFraction(min(exampleCalibrationTable$std_weights_kDa)),Inf)
  testthat::expect_gt(calibrationFunctions$MWtoFraction(min(exampleCalibrationTable$std_weights_kDa)),0)
  testthat::expect_lt(calibrationFunctions$FractionToMW(max(exampleCalibrationTable$std_elu_fractions)),Inf)
  testthat::expect_gt(calibrationFunctions$FractionToMW(min(exampleCalibrationTable$std_elu_fractions)),0)
  
})
