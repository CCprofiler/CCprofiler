context("generateHypotheses")


## Load example Data
complexHypotheses <- exampleComplexHypotheses

## Generate the binary network
binaryInteractions <- generateBinaryNetwork(complex_hypotheses = complexHypotheses)

## Calculate the path lengths
shortestPaths <- calculatePathlength(binaryInteractions)

## Generate a set of target hypotheses from the binary interactions
targetHypotheses <- generateComplexTargets(dist_info = shortestPaths,
                                            max_distance = 1,
                                            redundancy_cutoff = 1)
## Collapse redundancies
complexHypothesesCollapsed <- collapseHypothesis(targetHypotheses)

## Generate the Decoys
decoys <- generateComplexDecoys(target_hypotheses = complexHypotheses,
                                dist_info = shortestPaths,
                                min_distance = 2)


test_that("generateBinaryNetwork",{
  testthat::expect_error(generateBinaryNetwork(c()), "x is not a data.table")
  testthat::expect_identical(names(binaryInteractions), c("a", "b"))
  testthat::expect_false(any(duplicated(binaryInteractions)))
  #test if the interactions are number of subunits choose 2
  n_sub <- 5
  one_comp <- generateBinaryNetwork(data.table(complex_id = rep("1",n_sub), protein_id = 1:n_sub))
  testthat::expect_equal(nrow(one_comp), choose(n_sub, 2))
})

test_that("calculatePathlength",{
  testthat::expect_error(calculatePathlength(c()), "the data frame should contain at least two columns")
  testthat::expect_identical(names(shortestPaths), c("x", "y", "dist"))
  testthat::expect_false(any(duplicated(shortestPaths[,.(x,y)])))
  ## Test that no binary interactions get a distance of 2 or more
  testthat::expect_false(any(duplicated(rbind(shortestPaths[dist >= 2,.(x,y)], binaryInteractions[,.(x=a, y=b)]))))
})

test_that("generateComplexTargets",{
  testthat::expect_error(generateComplexTargets(c()), "non-character argument")
  testthat::expect_identical(names(targetHypotheses), c("complex_id", "complex_name", "protein_id"))
  testthat::expect_false(any(duplicated(targetHypotheses)))
  
  # Test that a single complex is made right
  one_prot <- shortestPaths[1,x]
  one_complex <- c(one_prot, shortestPaths[x == one_prot & dist == 1, y])
  targets <- copy(targetHypotheses)
  testthat::expect_true(any(targets[, all(protein_id %in% one_complex), by = complex_id]$V1))
})

test_that("collapseHypothesis",{
  testthat::expect_error(collapseHypothesis(c()))
  testthat::expect_identical(names(complexHypothesesCollapsed), c("complex_id", "complex_name", "protein_id"))
  testthat::expect_false(any(duplicated(complexHypothesesCollapsed)))
  ## Test specific scenarios
  n_sub <- 5
  comp1 <- data.table(complex_id = rep("1",n_sub), protein_id = as.character(1:n_sub))
  comp2 <- data.table(complex_id = rep("2",n_sub), protein_id = as.character(1:n_sub))
  testthat::expect_identical(collapseHypothesis(rbind(comp1, comp2))[,.(complex_id, protein_id)], comp1)
  testthat::expect_identical(collapseHypothesis(rbind(comp1, comp2), redundancy_cutoff = 0)[,.(complex_id, protein_id)], comp1)
  testthat::expect_identical(collapseHypothesis(rbind(comp1, comp2[1:3]))[,.(complex_id, protein_id)], comp1)
  testthat::expect_identical(collapseHypothesis(rbind(comp1, comp2[1:3]), redundancy_cutoff = 2)[,.(complex_id, protein_id)], comp1)
  
  comp3 <- data.table(complex_id = rep("2",n_sub), protein_id = as.character(4: (n_sub + 3)))
  testthat::expect_identical(collapseHypothesis(rbind(comp1, comp3), redundancy_cutoff = 0)[,.(complex_id, protein_id)], rbind(comp1, comp3))
  testthat::expect_identical(collapseHypothesis(rbind(comp1, comp3), redundancy_cutoff = 1)[,.(complex_id, protein_id)], rbind(comp1, comp3))
  testthat::expect_identical(collapseHypothesis(rbind(comp1, comp3), redundancy_cutoff = 1.3)[,.(complex_id, protein_id)], rbind(comp1, comp3))
  testthat::expect_identical(collapseHypothesis(rbind(comp1, comp3), redundancy_cutoff = 1.4)[,.(complex_id, protein_id)], 
                             rbind(comp1, data.table(complex_id = rep("1",3), protein_id = as.character(6:8))))

})

test_that("generateDecoys",{
  testthat::expect_error(generateDecoys(c()))
  testthat::expect_error(generateDecoys(c(),c(), c()))
  
  # Generate one decoy
  decoy <- generateDecoys(size = 3, all_proteins = unique(complexHypotheses$protein_id), 
                 dist_info = subset(shortestPaths, dist <= 2),
                 min_distance = 2,
                 n_tries = 3)
  decoy <-strsplit(decoy, ";")[[1]]
  testthat::expect_equal(length(decoy), 3)
  testthat::expect_true(all(shortestPaths[x %in% decoy & y %in% decoy & x != y, dist] > 2))
  
})

test_that("generateComplexDecoys",{
  testthat::expect_error(generateComplexDecoys(c()))
  testthat::expect_identical(names(decoys), c("complex_id", "complex_name", "protein_id"))
  testthat::expect_false(any(duplicated(decoys)))
  # Test size distribution
  testthat::expect_equal(nrow(complexHypotheses), nrow(decoys))
  testthat::expect_equal(sum(complexHypotheses[,.N, complex_id]$N), sum(decoys[,.N, complex_id]$N))
  testthat::expect_equal(sd(complexHypotheses[,.N, complex_id]$N), sd(decoys[,.N, complex_id]$N))
  # Test if parallelization has impact & seed
  testthat::expect_identical(decoys,
                             generateComplexDecoys(target_hypotheses = complexHypotheses,
                                                   dist_info = shortestPaths,
                                                   min_distance = 2, n_cores = 2)
  )
  testthat::expect_false(identical(decoys,
                                   generateComplexDecoys(target_hypotheses = complexHypotheses,
                                                         dist_info = shortestPaths,
                                                         min_distance = 2, seed = 234)
  ))
  #Test append
  testthat::expect_identical(rbind(complexHypotheses, decoys),
                             generateComplexDecoys(target_hypotheses = complexHypotheses,
                                                   dist_info = shortestPaths,
                                                   min_distance = 2, append = T)
  )
})


