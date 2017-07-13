context("featureMethods")

test_that("getBestFeatures",{
  ## data for tests
  protFeatures <- exampleProteinFeatures
  complexFeatures <- exampleComplexFeatures
  bestProteinFeatures <- getBestFeatures(protFeatures)
  bestComplexFeatures <- getBestFeatures(complexFeatures)
  ## test if only one feature per hypothesis
  testthat::expect_identical(length(unique(protFeatures$protein_id)), nrow(bestProteinFeatures))
  testthat::expect_identical(length(unique(complexFeatures$complex_id)), nrow(bestComplexFeatures))
  ## test if feature with highest number of subunits was selected as best feature
  testProteins <- unique(protFeatures$protein_id)[1:10]
  for (prot in testProteins) {
    testthat::expect_identical(max(subset(protFeatures,protein_id==prot)$n_subunits_detected), subset(bestProteinFeatures,protein_id==prot)$n_subunits_detected)
  }
  testComplexes <- unique(complexFeatures$protein_id)[1:10]
  for (complex in testComplexes) {
    testthat::expect_identical(max(subset(complexFeatures,complex_id==complex)$n_subunits_detected), subset(bestComplexFeatures,complex_id==complex)$n_subunits_detected)
  }
})

test_that("filterFeatures",{
  ## data for tests
  complexFeatures <- exampleComplexFeatures
  filteredComplexFeatures <- filterFeatures(feature_table=complexFeatures,
                                            complex_ids=NULL,
                                            protein_ids=NULL,
                                            min_feature_completeness=0.3,
                                            min_hypothesis_completeness=NULL,
                                            min_subunits=3,
                                            min_peak_corr=0.3,
                                            min_monomer_distance_factor=2)
  ## test completeness, n_subunits and peak_corr @TODO ids & min_monomer_distance
  testthat::expect_gte(min(filteredComplexFeatures$completeness),0.3)
  testthat::expect_gte(min(filteredComplexFeatures$n_subunits_detected),3)
  testthat::expect_gte(min(filteredComplexFeatures$peak_corr),0.3)
})

test_that("filterByStepwiseCompleteness",{
  ## data for tests
  complexFeatures <- exampleComplexFeatures
  filteredComplexFeatures <- filterByStepwiseCompleteness(feature_table=complexFeatures,
                                                          min_subunits_annotated=5,
                                                          completeness_vector=c(0.8,0.3),
                                                          level="hypothesis")
  filteredComplexFeaturesFeatures <- filterByStepwiseCompleteness(feature_table=complexFeatures,
                                                          min_subunits_annotated=5,
                                                          completeness_vector=c(0.8,0.3),
                                                          level="feature")
  ## Test data
  ## test hypothesis level
  sub_small <- subset(filteredComplexFeatures,n_subunits_annotated <= 5)
  sub_small <- getBestFeatures(sub_small)
  sub_big <- subset(filteredComplexFeatures,n_subunits_annotated > 5)
  sub_big<- getBestFeatures(sub_big)
  testthat::expect_gte(min(sub_small$completeness),0.8)
  testthat::expect_gte(min(sub_big$completeness),0.3)
  ## Test feature level
  sub_small_feature <- subset(filteredComplexFeaturesFeatures,n_subunits_annotated <= 5)
  sub_big_feature <- subset(filteredComplexFeaturesFeatures,n_subunits_annotated > 5)
  testthat::expect_gte(min(sub_small_feature$completeness),0.8)
  testthat::expect_gte(min(sub_big_feature$completeness),0.3)
})