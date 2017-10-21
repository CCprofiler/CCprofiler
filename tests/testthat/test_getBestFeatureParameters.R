context("getBestFeatureParameters")

#------------------------
## Complex level
#------------------------

## Load example data into lis to simulate grid search results
complexFeaturesGrid <- list(exampleComplexFeatures)

## Perform the filter grid search
complexFeaturesGridFiltered <- filterGridSearchResults(complexFeaturesGrid,
                                                       peak_corr_cutoffs = c(0.75,0.9),
                                                       feature_completeness_cutoffs = c(0.5,1),
                                                       hypothesis_completeness_cutoffs = c(0.5,1),
                                                       n_subunits_cutoffs =c(3,4),
                                                       monomer_distance_cutoffs = c(1,2),
                                                       remove_decoys=FALSE)
singleFilterC <- filterFeatures(exampleComplexFeatures,
                                min_feature_completeness = complexFeaturesGridFiltered[[1]]$min_feature_completeness[1], 
                                min_hypothesis_completeness = complexFeaturesGridFiltered[[1]]$min_hypothesis_completeness[1], 
                                min_subunits = complexFeaturesGridFiltered[[1]]$min_subunits[1],
                                min_peak_corr = complexFeaturesGridFiltered[[1]]$min_peak_corr[1], 
                                min_monomer_distance_factor = complexFeaturesGridFiltered[[1]]$min_monomer_distance_factor[1])

#------------------------
## Protein level
#------------------------

## Load example data into lis to simulate grid search results
proteinFeaturesGrid <- list(exampleProteinFeatures)

## Perform the filter grid search
proteinFeaturesGridFiltered <- filterGridSearchResults(proteinFeaturesGrid,
                                                       peak_corr_cutoffs = c(0.75,0.9),
                                                       feature_completeness_cutoffs = c(0.5,1),
                                                       hypothesis_completeness_cutoffs = c(0.5,1),
                                                       n_subunits_cutoffs =c(3,4),
                                                       monomer_distance_cutoffs = c(0,1),
                                                       remove_decoys=FALSE)
singleFilterP <- filterFeatures(exampleProteinFeatures,
                               min_feature_completeness = proteinFeaturesGridFiltered[[1]]$min_feature_completeness[1], 
                               min_hypothesis_completeness = proteinFeaturesGridFiltered[[1]]$min_hypothesis_completeness[1], 
                               min_subunits = proteinFeaturesGridFiltered[[1]]$min_subunits[1],
                               min_peak_corr = proteinFeaturesGridFiltered[[1]]$min_peak_corr[1], 
                               min_monomer_distance_factor = proteinFeaturesGridFiltered[[1]]$min_monomer_distance_factor[1])

test_that("filterGridSearchResults",{
  
  testthat::expect_identical(filterGridSearchResults(c()), NULL)
  
  
  testthat::expect_equal(length(complexFeaturesGridFiltered), 32)
  for(i in 1:length(complexFeaturesGridFiltered)){
    testthat::expect_identical(names(complexFeaturesGridFiltered[[i]]), 
                                     names(complexFeaturesGridFiltered[[1]]))
  } 
  colnamesC <- names(singleFilterP)
  testthat::expect_identical(proteinFeaturesGridFiltered[[1]][,colnamesC, with = F],
                             singleFilterP[,colnamesC, with = F])
  
  
  
  testthat::expect_equal(length(proteinFeaturesGridFiltered), 32)
  for(i in 1:length(proteinFeaturesGridFiltered)){
    testthat::expect_identical(names(proteinFeaturesGridFiltered[[i]]), 
                               names(proteinFeaturesGridFiltered[[1]]))
  }
  colnamesP <- names(singleFilterP)
  testthat::expect_identical(proteinFeaturesGridFiltered[[1]][,colnamesP, with = F],
                             singleFilterP[,colnamesP, with = F])
})

test_that("estimateGridSearchDecoyFDR",{
  
  grid_search_params =c("min_feature_completeness",
                        "min_hypothesis_completeness",
                        "min_subunits",
                        "min_peak_corr",
                        "min_monomer_distance_factor")
  gridStatsC <- estimateGridSearchDecoyFDR(complexFeaturesGridFiltered,
                                           grid_search_params)
  gridStatsP <- estimateGridSearchDecoyFDR(proteinFeaturesGridFiltered,
                                           grid_search_params)
  testthat::expect_error(estimateGridSearchDecoyFDR(c()), "subscript out of bounds")
  
  testthat::expect_equal(ncol(gridStatsC), length(grid_search_params) + 4)
  testthat::expect_equal(ncol(gridStatsP), length(grid_search_params) + 4)
  testthat::expect_equal(as.numeric(gridStatsC[1,1, drop = T]),
                         estimateDecoyFDR(complexFeaturesGridFiltered[[1]])[[1]])
  testthat::expect_equal(as.numeric(gridStatsP[1,1, drop = T]),
                         estimateDecoyFDR(proteinFeaturesGridFiltered[[1]])[[1]])
  
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

