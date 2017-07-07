context("importPCPdata")



## Sum the intensities of the top 2 peptides to get protein intensities
protTracesFilt <- proteinQuantification(examplePeptideTracesFiltered,
                                    topN = 2)

protTraces <- proteinQuantification(examplePeptideTraces,
                                    topN = 2)
examplePeptideTracesSpc <- calculateSibPepCorr(examplePeptideTraces,
                                               plot = F,
                                               PDF = F)
protTracesSpc <- proteinQuantification(examplePeptideTracesSpc,
                                    topN = 2)

test_that("Fringe cases", {
  testthat::expect_error(proteinQuantification(c()),"Object is not of class traces.")
})

test_that("Output format",{
  ## Test that sibPepCorr is handelled correctly
  testthat::expect_identical(protTracesSpc$traces, protTraces$traces)
  testthat::expect_identical(protTracesSpc$fraction_annotation, protTraces$fraction_annotation)
  testthat::expect_identical(protTracesSpc$trace_type, protTraces$trace_type)
  testthat::expect_identical(protTracesSpc$trace_annotation[, !"SibPepCorr_protein_mean", with = FALSE],
                             protTraces$trace_annotation)
  testthat::expect_equal(examplePeptideTracesSpc$trace_annotation[protein_id == "Q9Y316", mean(SibPepCorr)],
                         protTracesSpc$trace_annotation[protein_id == "Q9Y316", SibPepCorr_protein_mean])
  ## Test against example data
  testthat::expect_equal(protTracesFilt, exampleProteinTraces)
})
