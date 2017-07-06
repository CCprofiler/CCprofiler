context("tracesMethods")

test_that("subset.traces",{
  ## get data for tests
  subsetPeptides <- c("AIIDEFEQK","AIQLSGAEQLEALK","AKEALIAASETLK")
  subsetProtein <- "Q15021"
  subsetFractions <- c(5:20)
  peptideSubsettedPeptideTraces <- subset(examplePeptideTraces,trace_subset_ids=subsetPeptides,fraction_ids=subsetFractions)
  proteinSubsettedPeptideTraces <- subset(examplePeptideTraces,trace_subset_ids=subsetProtein,trace_subset_type="protein_id")
  proteinSubsettedProteinTraces <- subset(exampleProteinTraces,trace_subset_ids=subsetProtein)
  ## test if subsetting produced valid traces object.
  testthat::expect_null(.tracesTest(peptideSubsettedPeptideTraces))
  testthat::expect_null(.tracesTest(proteinSubsettedPeptideTraces))
  testthat::expect_null(.tracesTest(proteinSubsettedProteinTraces))
  ## test if peptide id subsetting worked fine
  testthat::expect_equal(length(peptideSubsettedPeptideTraces$traces$id), length(subsetPeptides))
  testthat::expect_equal(length(peptideSubsettedPeptideTraces$trace_annotation$id), length(subsetPeptides))
  ## test if fraction subsetting worked fine
  testthat::expect_equal(length(names(peptideSubsettedPeptideTraces$traces)), length(subsetFractions)+1)
  testthat::expect_equal(peptideSubsettedPeptideTraces$fraction_annotation$id, subsetFractions)
  ## test proteinTraces subsetting
  testthat::expect_equal(length(proteinSubsettedProteinTraces$traces$id), length(subsetProtein))
  testthat::expect_equal(length(proteinSubsettedProteinTraces$trace_annotation$id), length(subsetProtein))
  ## test if subsetting by thomething other than "id" works
  testthat::expect_equal(length(proteinSubsettedPeptideTraces$traces$id), 42)
  testthat::expect_equal(length(proteinSubsettedPeptideTraces$trace_annotation$id), 42)
  ## test if traces and trace_annotation are identical (no shift in order etc)
  testthat::expect_identical(proteinSubsettedPeptideTraces$trace_annotation$id,proteinSubsettedPeptideTraces$traces$id)
  })

test_that("toLongFormat",{
  ## get data for tests
  tracesWide <- examplePeptideTraces$traces
  tracesLong <- toLongFormat(tracesWide)
  ## test
  testthat::expect_identical(unique(tracesWide$id),unique(tracesLong$id))
  testthat::expect_identical(c(unique(tracesLong$fraction),"id"),names(tracesWide))
  })

test_that("annotateMolecularWeight",{
  calibration <- calibrateSECMW(exampleCalibrationTable)
  mwTraces <- annotateMolecularWeight(examplePeptideTraces, calibration)
  testthat::expect_null(.tracesTest(mwTraces))
  testthat::expect_true("molecular_weight" %in% names(mwTraces$fraction_annotation))
  })

test_that("summary.traces",{
  tracesSummary = summary(examplePeptideTraces)
  testthat::expect_identical(tracesSummary$type,"peptide")
  testthat::expect_equal(tracesSummary$fraction_count,81)
  testthat::expect_identical(tracesSummary$annotations,
    c("protein_id","id","Entry_name","Status","Protein_names","Gene_names","Organism","Length","Mass","protein_mw"))
  testthat::expect_equal(as.numeric(tracesSummary$metrics),c(1743,1506,237,14))
  })
