context("Import from OpenSWATH")


traces <- importFromOpenSWATH(data = exampleOpenSWATHinput,
                              annotation_table = exampleFractionAnnotation,
                              rm_requantified = TRUE,verbose=FALSE)
tracesFromWide <- importPCPdata(input_data = examplePCPdataWide,
                                fraction_annotation = exampleFractionAnnotation)

test_that("Data import messages", {
  testthat::expect_error(importFromOpenSWATH(),"Need to specify data in form of OpenSWATH result file or R data.table.")
  testthat::expect_error(importFromOpenSWATH(data=exampleOpenSWATHinput),"Need to specify annotation_table.")
  testthat::expect_error(importFromOpenSWATH(data="mock",annotation_table=exampleFractionAnnotation),"data file doesn't exist")
  mock <- c()
  testthat::expect_error(importFromOpenSWATH(data=mock,annotation_table=exampleFractionAnnotation),"data input is neither file name or data.table")
  testthat::expect_error(importFromOpenSWATH(data=exampleOpenSWATHinput,annotation_table="mock"),"annotation_table file doesn't exist")
  testthat::expect_error(importFromOpenSWATH(data=exampleOpenSWATHinput,annotation_table=head(exampleFractionAnnotation)),"Number of file names in annotation_table does not match data")
  testthat::expect_error(importFromOpenSWATH(data=head(exampleOpenSWATHinput),annotation_table=exampleFractionAnnotation),"Number of file names in annotation_table does not match data")
})

test_that("Output format",{
  testthat::expect_equal(sum(subset(traces$traces, select = -id)), sum(subset(tracesFromWide$traces, select = -id)))
})
