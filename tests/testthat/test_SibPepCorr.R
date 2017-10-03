context("Sbibling Peptide Correlation Functions")

FFT <- 0.4
tracesRawSpc <- calculateSibPepCorr(traces = examplePeptideTraces,
                                    plot = FALSE)
                                    
tracesFiltered <- filterBySibPepCorr(traces = tracesRawSpc,
                                     fdr_cutoff = 0.01,
                                     fdr_type = "protein",
                                     FFT = FFT,
                                     absolute_spcCutoff = NULL,
                                     rm_decoys = FALSE,
                                     plot = FALSE,
                                     CSV = FALSE)

tracesFilteredPep <- filterBySibPepCorr(traces = tracesRawSpc,
                                     fdr_cutoff = 0.01,
                                     fdr_type = "peptide",
                                     FFT = FFT,
                                     absolute_spcCutoff = NULL,
                                     rm_decoys = FALSE,
                                     plot = FALSE,
                                     CSV = FALSE)

FDRtable <- rocSibPepCorr(traces = tracesRawSpc,
                          plot = FALSE, 
                          fdr_type = "peptide",
                          stepsize = 0.01,
                          FFT = FFT)

spc_cutoff <- FDRtable[50]$SibPepCorr_cutoff

test_that("Fringe cases", {
  testthat::expect_error(filterBySibPepCorr(c()),"Object is not of class traces.")
  testthat::expect_error(filterBySibPepCorr(exampleProteinTraces),"Traces object is of wrong type. Please check your input traces.")
  testthat::expect_error(importPCPdata(rbind(examplePCPdataWide[1], examplePCPdataWide[1]), exampleFractionAnnotation))
  testthat::expect_error(importPCPdata(rbind(exampleDataProtWide[1], exampleDataProtWide[1]), exampleFractionAnnotation))
  
})

test_that("Output format",{
  ## Test if ountouched objects are the same
  testthat::expect_identical(tracesFiltered$trace_type, examplePeptideTraces$trace_type)
  testthat::expect_identical(tracesFiltered$fraction_annotation, examplePeptideTraces$fraction_annotation)
  ## Test if FDR calculations are sensical
  testthat::expect_equal(FDRtable[1,n_targetPeptides] + FDRtable[1,n_decoyPeptides], 
                         nrow(tracesRawSpc$trace_annotation[!is.na(SibPepCorr)]))
  testthat::expect_true(all(FDRtable$proteinFDR >= 0 & FDRtable$proteinFDR <= 1))
  testthat::expect_true(all(FDRtable$peptideFDR >= 0 & FDRtable$peptideFDR <= 1))
  testthat::expect_equal(FDRtable$proteinFDR, FFT * FDRtable$n_decoyProteins / FDRtable$n_targetProteins)
  testthat::expect_equal(FDRtable$peptideFDR, FFT * FDRtable$n_decoyPeptides / FDRtable$n_targetPeptides)
  testthat::expect_equal(FDRtable$n_true_targetProteins, as.integer(FDRtable$n_targetProteins * (1-FDRtable$proteinFDR)))
  testthat::expect_equal(FDRtable$n_true_targetPeptides, as.integer(FDRtable$n_targetPeptides * (1-FDRtable$peptideFDR)))
  testthat::expect_lt(min(tracesFilteredPep$trace_annotation$SibPepCorr),
                         FDRtable[(peptideFDR <= 0.01), min(SibPepCorr_cutoff)])
})
