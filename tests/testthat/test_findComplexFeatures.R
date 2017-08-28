context("findComplexFeatures")

test_that("imputeMissingValues",{
  tracesMatrix <- getIntensityMatrix(exampleProteinTraces)
  tracesImputed_5p <- imputeMissingValues(tracesMatrix,"5%")
  tracesImputed_val <- imputeMissingValues(tracesMatrix,1587)
  testthat::expect_equal(tracesImputed_5p,tracesImputed_val)
  testthat::expect_equal(sum(tracesImputed_5p==0),0)
})


test_that("collapseComplexFeatures",{
  testDataCollapsing <- data.table(
    subgroup = c("Q15021;Q15003","Q15021;Q15003",
                 "Q9BPX3;Q15021","Q9BPX3;Q15021;Q15003",
                 "Q9NTJ3;O95347","Q9NTJ3;Q15003",
                 "Q9NTJ3;Q15003;O95347","Q9NTJ3;Q15021","Q9NTJ3;Q15021;Q15003"),
    left_sw = c(8,66,3,1,1,13,9,20,16),
    right_sw = c(27,81,18,17,26,30,25,36,34),
    score = c(0.9576664,0.9569117,0.9573860,0.9835337,0.9954245,0.9678766,0.9757193,0.9577040,0.9684879),
    n_subunits = c(2,2,2,3,2,2,3,2,3),
    apex = c(10,NA,10,10,9,NA,9,NA,NA),
    left_pp = c(1,NA,1,1,1,NA,1,NA,NA),
    right_pp = c(16,NA,16,16,16,NA,16,NA,NA),
    area = NA
  )
  
  testCollapse5 <- data.table(
    subgroup = c("Q15021;Q15003;Q9BPX3;Q9NTJ3;O95347"),
    left_sw = c(1),
    right_sw = c(17),
    score = c(0.9835337),
    n_subunits = c(5),
    apex = c(10),
    left_pp = c(1),
    right_pp = c(16),
    area = NA
  )
  
  testCollapse0 <- data.table(
    subgroup = c("Q15021;Q15003;Q9BPX3","Q9NTJ3;O95347;Q15003"),
    left_sw = c(1,9),
    right_sw = c(17,25),
    score = c(0.9835337,0.9757193),
    n_subunits = c(3,3),
    apex = c(10,9),
    left_pp = c(1,1),
    right_pp = c(16,16),
    area = NA
  )
  
  collapse5_apex <- collapseComplexFeatures(testDataCollapsing,rt_height=5,collapse_method="apex_only")
  collapse0_apex <- collapseComplexFeatures(testDataCollapsing,rt_height=0,collapse_method="apex_only")
  collapse5_network <- collapseComplexFeatures(testDataCollapsing,rt_height=5,collapse_method="apex_network")
  collapse0_network <- collapseComplexFeatures(testDataCollapsing,rt_height=0,collapse_method="apex_network")
  
  testthat::expect_equal(collapse5_apex,collapse5_network)
  testthat::expect_equal(collapse0_apex,collapse0_network)
  testthat::expect_equal(collapse5_apex,testCollapse5)
  testthat::expect_equal(collapse0_apex,testCollapse0)
  
  testDataCollapsing_diff <- data.table(
    subgroup = c("Q13330;Q09028","Q14839;Q09028",
                 "Q14839;Q13330;Q09028","Q14839;Q13330;Q09028",
                 "Q92769;Q13547","Q92769;Q13547","Q92769;Q13547;O60341"),
    left_sw = c(12,1,2,2,1,1,16),
    right_sw = c(39,26,21,21,39,39,31),
    score = c(0.9522825,0.9815125,0.9750086,0.9750086,0.9833626,0.9833626,0.9613223),
    n_subunits = c(2,2,3,3,2,2,3),
    apex = c(17,17,17,7,17,9,17),
    left_pp = c(11,11,11,1,12,1,14),
    right_pp = c(22,22,22,11,22,12,22),
    area = NA
  )
  
  test_collapse5_diff_apex <- data.table(
    subgroup = c("Q13330;Q09028;Q14839;Q92769;Q13547;O60341",
                 "Q14839;Q13330;Q09028;Q92769;Q13547"),
    left_sw = c(2,2),
    right_sw = c(21,21),
    score = c(0.9750086,0.9750086),
    n_subunits = c(6,5),
    apex = c(17,7),
    left_pp = c(11,1),
    right_pp = c(22,11),
    area = NA
  )
  
  test_collapse0_diff_apex <- data.table(
    subgroup = c("Q13330;Q09028;Q14839;Q92769;Q13547;O60341",
                 "Q14839;Q13330;Q09028","Q92769;Q13547"),
    left_sw = c(2,2,1),
    right_sw = c(21,21,39),
    score = c(0.9750086,0.9750086,0.9833626),
    n_subunits = c(6,3,2),
    apex = c(17,7,9),
    left_pp = c(11,1,1),
    right_pp = c(22,11,12),
    area = NA
  )
  
  test_collapse5_diff_network <- data.table(
    subgroup = c("Q13330;Q09028;Q14839","Q92769;Q13547;O60341","Q14839;Q13330;Q09028","Q92769;Q13547"),
    left_sw = c(2,16,2,1),
    right_sw = c(21,31,21,39),
    score = c(0.9750086,0.9613223,0.9750086,0.9833626),
    n_subunits = c(3,3,3,2),
    apex = c(17,17,7,9),
    left_pp = c(11,14,1,1),
    right_pp = c(22,22,11,12),
    area = NA
  )
  
  collapse5_diff_apex <- collapseComplexFeatures(testDataCollapsing_diff,rt_height=5,collapse_method="apex_only")
  collapse0_diff_apex <- collapseComplexFeatures(testDataCollapsing_diff,rt_height=0,collapse_method="apex_only")
  collapse5_diff_network <- collapseComplexFeatures(testDataCollapsing_diff,rt_height=5,collapse_method="apex_network")
  collapse0_diff_network <- collapseComplexFeatures(testDataCollapsing_diff,rt_height=0,collapse_method="apex_network")
  
  testthat::expect_equal(collapse5_diff_network,collapse0_diff_network)
  testthat::expect_equal(collapse5_diff_apex,test_collapse5_diff_apex)
  testthat::expect_equal(collapse0_diff_apex,test_collapse0_diff_apex)
  testthat::expect_equal(collapse5_diff_network,test_collapse5_diff_network)
})

test_that("calculateFeatureCorrelation",{
  testData <- data.table(
    subgroup = c("Q13330;Q09028;Q14839","Q92769;Q13547;O60341","Q14839;Q13330;Q09028","Q92769;Q13547"),
    left_sw = c(2,16,2,1),
    right_sw = c(21,31,21,39),
    score = c(0.9750086,0.9613223,0.9750086,0.9833626),
    n_subunits = c(3,3,3,2),
    apex = c(17,17,7,9),
    left_pp = c(11,14,1,1),
    right_pp = c(22,22,11,12),
    area = NA
  )
  tracesMatrix <- getIntensityMatrix(exampleProteinTraces)
  imputedMatrix <- imputeMissingValues(tracesMatrix,"5%")
  testDataCorr <- calculateFeatureCorrelation(imputedMatrix, testData, toTable = FALSE)
  testthat::expect_equal(testData[,c(1:8)],testDataCorr[,c(1:8)])
  testthat::expect_equal(testDataCorr$area,c(1415864,808581,409433,499226))
  testthat::expect_equal(testDataCorr$peak_corr,c(0.9547027,0.9573859,0.7322601,0.9796383), tolerance = .0000001)
})

test_that("estimateComplexFeatureStoichiometry",{
  testData <- data.table(
    subgroup = c("Q13330;Q09028;Q14839","Q92769;Q13547;O60341","Q14839;Q13330;Q09028","Q92769;Q13547","Q92769;Q13547"),
    left_sw = c(2,16,2,1,1),
    right_sw = c(21,31,21,39,39),
    score = c(0.9750086,0.9613223,0.9750086,0.9833626,0.9833626),
    n_subunits = c(3,3,3,2,2),
    apex = c(17,17,7,9,9),
    left_pp = c(11,14,1,1,1),
    right_pp = c(22,22,11,12,12),
    area = NA
  )
  tracesMatrix <- getIntensityMatrix(exampleProteinTraces)
  imputedMatrix <- imputeMissingValues(tracesMatrix,"5%")
  testDataCorr <- calculateFeatureCorrelation(imputedMatrix, testData, toTable = FALSE)
  testStoichiometry <- data.table(
    left_sw = c(2,16,2,1,1),
    right_sw = c(21,31,21,39,39),
    score = c(0.9750086,0.9613223,0.9750086,0.9833626,0.9833626),
    n_subunits = c(3,3,3,2,2),
    apex = c(17,17,7,9,9),
    left_pp = c(11,14,1,1,1),
    right_pp = c(22,22,11,12,12),
    area = c(1415864,808581,409433,499226,499226),
    peak_corr = c(0.9547027,0.9573859,0.7322601,0.9796383,0.9796383),
    id = c("Q09028;Q13330;Q14839","O60341;Q13547;Q92769","Q09028;Q13330;Q14839","Q13547;Q92769","Q13547;Q92769"),
    total_intensity = c("695667;369472;350725","242329;297031;269221","234247;81089;94030","276525;222701","276525;222701"),
    intensity_ratio = c("1.98351129802552;1.0534521348635;1",
                        "1;1.22573443541631;1.1109730985561",
                        "2.88876419736339;1;1.15959008003552",
                        "1.24168728474502;1","1.24168728474502;1"),
    stoichiometry=c("2;1;1","1;1;1","3;1;1","1;1","1;1")
  )
  stoichiometryData <- estimateComplexFeatureStoichiometry(exampleProteinTraces,testDataCorr)
  testthat::expect_equal(testData[,c(2:8)],stoichiometryData[,c(1:7)])
  testthat::expect_equal(stoichiometryData,testStoichiometry, tolerance = .0000001)
})

test_that("annotateComplexFeatures",{
  testData <- data.table(
    left_sw = c(2,16,2,1,1),
    right_sw = c(21,31,21,39,39),
    score = c(0.9750086,0.9613223,0.9750086,0.9833626,0.9833626),
    n_subunits = c(3,3,3,2,2),
    apex = c(17,17,7,9,9),
    left_pp = c(11,14,1,1,1),
    right_pp = c(22,22,11,12,12),
    area = c(1415864,808581,409433,499226,499226),
    peak_corr = c(0.9547027,0.9573859,0.7322601,0.9796383,0.9796383),
    id = c("Q09028;Q13330;Q14839","O60341;Q13547;Q92769","Q09028;Q13330;Q14839","Q13547;Q92769","Q13547;Q92769"),
    total_intensity = c("695667;369472;350725","242329;297031;269221","234247;81089;94030","276525;222701","276525;222701"),
    intensity_ratio = c("1.98351129802552;1.0534521348635;1",
                        "1;1.22573443541631;1.1109730985561",
                        "2.88876419736339;1;1.15959008003552",
                        "1.24168728474502;1","1.24168728474502;1"),
    stoichiometry=c("2;1;1","1;1;1","3;1;1","1;1","1;1")
  )  
  complex.subunits <- exampleComplexHypotheses[complex_id == "614",
                                              protein_id]
  traces.subs <- subset(traces=exampleProteinTraces,trace_subset_ids=complex.subunits)
  complex.annotation <- subset(exampleComplexHypotheses,complex_id == "614")
  
  annotatedData <- annotateComplexFeatures(traces.subs,testData,complex.annotation)
  testthat::expect_equal(nrow(testData),nrow(annotatedData))
})