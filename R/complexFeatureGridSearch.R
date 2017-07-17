#' Perform complex feature grid search
#' @description Perform complex feature finding (calls \code{\link{findComplexFeatures}})
#' with all possible combinations of the specified parameters (grid search).
#' @details The runtime of this function scales with the binomial coefficient of the total 
#' number of parameters specified, and can therefore take a long time. If many parameter combinations
#' are searched parallelization is strongly recommended.
#' @param traces traces object of type protein.
#' @param corrs Numeric, vector of correlation_cutoff values to test.
#' @param windows Numeric, vector of window_size values to test.
#' @param smoothing Numeric, vector of smoothing_length values to test.
#' @param rt_heights Numeric, vector of rt_height values to test.
#' @param n_cores Numeric, number of cores to use (default=1).
#' @return List of search result tables for every possible parameter combination.
#' The result tables contain additional columns specifying the parameters.
#' @export
#' @examples 
#' 
#' ## Load example data
#' proteinTraces <- exampleProteinTraces
#' complexHypotheses <- exampleComplexHypotheses
#' 
#' ## Perform a small grid search for 2 parameter combinations
#' gridList <- performComplexGridSearch(traces = proteinTraces,
#'                          complex_hypothesis = complexHypotheses,
#'                          corrs = c(0.5, 0.9),
#'                          windows = 10,
#'                          smoothing = 7,
#'                          rt_heights = 4,
#'                          n_cores = 2)
#'                          
#' lapply(gridList, head, n = 2)

performComplexGridSearch <- function(traces,
                                     complex_hypothesis,
                                     corrs = c(0.5,0.75,0.9,0.95),
                                     windows = c(8,10,12),
                                     smoothing = c(7,9,11),
                                     rt_heights = c(3,4,5),
                                     n_cores=1){
  
  .tracesTest(traces, type = "protein")
  .testGridParameter(corrs, "corrs")
  .testGridParameter(windows, "windows")
  .testGridParameter(smoothing, "smoothing")
  .testGridParameter(rt_heights, "rt_heights")
  
  parameter_grid <- as.data.table(expand.grid(corrs,windows,smoothing,rt_heights))
  names(parameter_grid) <- c("corr","window","smoothing","rt_height")
  if(n_cores > nrow(parameter_grid)){
    n_cores <- nrow(parameter_grid)
    message(paste0("Only ",n_cores," required for parallelization. Using ",n_cores," cores."))
  }
  
  cl <- snow::makeCluster(n_cores)
  # setting a seed is absolutely crutial to ensure reproducible results!
  clusterSetRNGStream(cl,123)
  doSNOW::registerDoSNOW(cl)
  clusterEvalQ(cl,library(SECprofiler,data.table))
  data <- parRapply(cl,parameter_grid,
                    FUN = .runGridComplexFeatureFinding,
                    protTraces = traces,
                    complex_hypothesis = complex_hypothesis)
  stopCluster(cl)
  
  return(data)
}

.runGridComplexFeatureFinding <- function(params,
                                          protTraces,
                                          complex_hypothesis) {
  res <- findComplexFeatures(traces=protTraces,
                             complex_hypothesis = complex_hypothesis,
                             corr_cutoff = as.numeric(params["corr"]),
                             window_size = as.numeric(params["window"]),
                             parallelized = FALSE,
                             collapse_method="apex_network",
                             perturb_cutoff="1%",
                             rt_height=as.numeric(params["rt_height"]),
                             smoothing_length=as.numeric(params["smoothing"]))
  
  res$corr <- as.numeric(params["corr"])
  res$window <- as.numeric(params["window"])
  res$rt_height <- as.numeric(params["rt_height"])
  res$smoothing_length <- as.numeric(params["smoothing"])
  #saveRDS(res,file=paste0("complexFeatures/complexFeatures_",name,".rda"))
  res[]
}

