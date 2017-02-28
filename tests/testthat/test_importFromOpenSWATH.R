context("Import from OpenSWATH")

#@TODO put test data in lazy load
test_data='/Volumes/ibludau-1/SEC/write_tests/swath_import_reduced.tsv'
test_annotation_table="/Volumes/ibludau-1/SEC/write_tests/annotation_SECprofiler_noInput.txt"
output_rds=readRDS("/Volumes/ibludau-1/SEC/write_tests/test_osw_input.rda")

test_that("Data import messages", {
  expect_error(importFromOpenSWATH(),"Need to specify data in form of OpenSWATH result file or R data.table.")
  expect_error(importFromOpenSWATH(data=test_data),"Need to specify annotation_table.")
  expect_error(importFromOpenSWATH(data="mock",annotation_table=test_annotation_table),"data file doesn't exist")
  mock <- c()
  expect_error(importFromOpenSWATH(data=mock,annotation_table=test_annotation_table),"data input is neither file name or data.table")
  expect_error(importFromOpenSWATH(data=test_data,annotation_table="mock"),"annotation_table doesn't exist")
  expect_error(importFromOpenSWATH(data=test_data,annotation_table=head(test_annotation_table)),"Number of file names in annotation does not match data")
  expect_error(importFromOpenSWATH(data=head(test_data),annotation_table=test_annotation_table),"Number of file names in annotation does not match data")
})

test_that("Output format",{
  expect_equal(importFromOpenSWATH(data=test_data,annotation_table=test_annotation_table),output_rds)
  decoy_idx=grep("DECOY",output_rds$trace_annotation$protein_id,invert=TRUE)
  output_rds_noDecoys=output_rds
  output_rds_noDecoys$trace_annotation=output_rds_noDecoys$trace_annotation[decoy_idx]
  output_rds_noDecoys$traces=output_rds_noDecoys$traces[decoy_idx]
  expect_equal(importFromOpenSWATH(data=test_data,annotation_table=test_annotation_table,rm_decoys=TRUE),output_rds_noDecoys)
  #@TODO write test for MS1Quant
  #@TODO write test for import of non-uniprot names
  #@TODO write test for rm_requantified
})
