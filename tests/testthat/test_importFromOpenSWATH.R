context("Import from OpenSWATH")

#@TODO put test data in lazy load

output_rds=readRDS("/Volumes/ibludau-1/SEC/write_tests/test_osw_input.rda")

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
  testthat::expect_equal(importFromOpenSWATH(data=exampleOpenSWATHinput,annotation_table=exampleFractionAnnotation),output_rds)
  decoy_idx=grep("DECOY",output_rds$trace_annotation$protein_id,invert=TRUE)
  output_rds_noDecoys=output_rds
  output_rds_noDecoys$trace_annotation=output_rds_noDecoys$trace_annotation[decoy_idx]
  output_rds_noDecoys$traces=output_rds_noDecoys$traces[decoy_idx]
  testthat::expect_equal(importFromOpenSWATH(data=exampleOpenSWATHinput,annotation_table=exampleFractionAnnotation,rm_decoys=TRUE),output_rds_noDecoys)
  #@TODO write test for MS1Quant
  #@TODO write test for import of non-uniprot names
  #@TODO write test for rm_requantified
})
