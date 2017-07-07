context("calibrateMW")

calibrationFunctions <- calibrateSECMW(exampleCalibrationTable,
                                       plot = T,
                                       PDF = FALSE)

test_that("Fringe cases", {
  testthat::expect_error(calibrateSECMW(c()),"calibration_table input is neither file name or data.table")
  testthat::expect_error(calibrateSECMW(data.table()),"calibration_table seems to be empty")
  
})
test_that("Output Format", {
  testthat::expect_equal(calibrationFunctions$MWtoSECfraction(0), Inf)
  testthat::expect_lt(min(exampleCalibrationTable$std_weights_kDa),Inf)
  testthat::expect_gt(min(exampleCalibrationTable$std_weights_kDa),0)
  testthat::expect_lt(max(exampleCalibrationTable$std_weights_kDa),Inf)
  testthat::expect_gt(max(exampleCalibrationTable$std_weights_kDa),0)
  
})
