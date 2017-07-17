
complexFeaturesUnique <- getUniqueFeatureGroups(exampleComplexFeatures,
                                                rt_height=2,
                                                distance_cutoff=1)

context("getUniqueFeatureGroups")
test_that("Output format",{
  testEqual <- copy(complexFeaturesUnique)
  testEqual <- testEqual[,unique_feature_identifier:=NULL]
  testthat::expect_true(all.equal(testEqual, exampleComplexFeatures, ignore.row.order=TRUE))
  testthat::expect_true("unique_feature_identifier" %in% names(complexFeaturesUnique))
})

context("callapseByUniqueFeatureGroups")
test_that("Output format",{
  complexFeaturesCollapsed <- callapseByUniqueFeatureGroups(complexFeaturesUnique,
                                                            rm_decoys=TRUE)
  testthat::expect_true(all(unique(complexFeaturesCollapsed$unique_feature_identifier) %in% unique(complexFeaturesUnique$unique_feature_identifier)))
  complexFeaturesCollapsed[ , `:=`( COUNT = .N , IDX = 1:.N ) , by = unique_feature_identifier ]
  testthat::expect_equal(max(complexFeaturesCollapsed$COUNT),1)
})