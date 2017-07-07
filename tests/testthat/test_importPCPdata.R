context("importPCPdata")



tracesFromLong <- importPCPdata(input_data = examplePCPdataLong,
                                fraction_annotation = exampleFractionAnnotation)
tracesFromWide <- importPCPdata(input_data = examplePCPdataWide,
                                fraction_annotation = exampleFractionAnnotation)
tracesFromWide_nodecoy <- importPCPdata(input_data = examplePCPdataWide,
                                        fraction_annotation = exampleFractionAnnotation,
                                        rm_decoys = TRUE)

exampleDataProt <- examplePCPdataLong[, .(intensity = sum(intensity)), by = .(protein_id, filename)]
exampleDataProtWide <- dcast(exampleDataProt, protein_id ~ filename, value.var = "intensity")
tracesFromLongProt <- importPCPdata(input_data = exampleDataProt,
                                    fraction_annotation = exampleFractionAnnotation)

examplePCPdataWide_scr <- cbind(examplePCPdataWide[,-15], examplePCPdataWide[,15])
tracesFromWide_scr <- importPCPdata(input_data = examplePCPdataWide_scr,
                                    fraction_annotation = exampleFractionAnnotation)

test_that("Fringe cases", {
  testthat::expect_error(importPCPdata(c()),"Need to specify fraction_annotation in form of a filename or R data.table.")
  testthat::expect_error(importPCPdata(c(), c()),"input_data input is neither filename or data.table")
  testthat::expect_error(importPCPdata(examplePCPdataLong, c()),"fraction_annotation input is neither filename or data.table")
  testthat::expect_error(importPCPdata(c(), exampleFractionAnnotation),"input_data input is neither filename or data.table")
  testthat::expect_error(importPCPdata(examplePCPdataLong, exampleFractionAnnotation[-1]),"Number of filenames in fraction_annotation does not match input_data.")
  testthat::expect_error(importPCPdata(rbind(examplePCPdataWide[1], examplePCPdataWide[1]), exampleFractionAnnotation))
  testthat::expect_error(importPCPdata(rbind(exampleDataProtWide[1], exampleDataProtWide[1]), exampleFractionAnnotation))
  
})

test_that("Output format",{
  testthat::expect_equal(tracesFromLong, tracesFromWide)
  # testthat::expect_equal(tracesFromLong, examplePeptideTracesUnannotated)
  testthat::expect_equal(tracesFromWide, tracesFromWide_scr)
  testthat::expect_equal(length(grep("DECOY", tracesFromWide_nodecoy$trace_annotation$protein_id)), 0)
  
})
