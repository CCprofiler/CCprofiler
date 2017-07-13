
#' Perform Protein feature grid search
#' @description Perform Protein feature finding (calls \code{\link{findProteinFeatures}})
#' with all possible combinations of the specified parameters (grid search).
#' @details The runtime of this function scales with the binomial coefficient of the total 
#' number of parameters specified, and can therefore take a long time. If many parameter combinations
#' are searched parallelization is strongly recommended.
#' @import data.table
#' @param traces traces object of type peptide
#' @param calibration list of two functions for calibration
#' @param corrs numeric vector
#' @param windows numeric vector
#' @param smoothing numeric vector
#' @param rt_heights numeric vector
#' @param parallelized logical default=TRUE
#' @param n_cores numeric number of cores to use if parallelized (default=1)
#' @return List with stats
#' @export
#' @examples  
#' 
#' ## Load example data and subset to reduce runtime
#' peptideTraces <- examplePeptideTraces
#' peptideTraces <- subset(peptideTraces,
#'                         trace_subset_ids = unique(peptideTraces$trace_annotation$protein_id)[2:3],
#'                         trace_subset_type = "protein_id")
#' 
#' ## Perform a small grid search for 2 parameter combinations
#' # Depending on the computational resources this can take several minutes 
#' gridList <- performProteinGridSearch(traces = peptideTraces,
#'                                      corrs = c(0.5, 0.9),
#'                                      windows = 10,
#'                                      smoothing = 7,
#'                                      rt_heights = 4,
#'                                      n_cores = 4)
#' 

performProteinGridSearch <- function(traces,
                                    corrs = c(0.5,0.75,0.9,0.95),
                                    windows = c(8,10,12),
                                    smoothing = c(7,9,11),
                                    rt_heights = c(3,4,5),
                                    parallelized = TRUE,
                                    n_cores=1
                                    ){
  parameter_grid <- expand.grid(corrs,windows,smoothing,rt_heights)
  names(parameter_grid) <- c("corr","window","smoothing","rt_height")
  cl <- snow::makeCluster(n_cores)
  # setting a seed is absolutely crutial to ensure reproducible results!
  clusterSetRNGStream(cl,123)
  doSNOW::registerDoSNOW(cl)
  clusterEvalQ(cl,library(SECprofiler,data.table))
  clusterExport(cl, list("generateRandomPepTraces"))
  data <- list()
  data <- c(data,parRapply(cl,parameter_grid,
                           FUN = .runGridProteinFeatureFinding,
                           traces = traces))
  stopCluster(cl)
  data
}

.runGridProteinFeatureFinding <- function(params,traces) {
  res = findProteinFeatures(traces=traces,
                            corr_cutoff = as.numeric(params["corr"]),
                            window_size = as.numeric(params["window"]),
                            parallelized = FALSE,
                            collapse_method="apex_only",
                            perturb_cutoff="1%",
                            rt_height=as.numeric(params["rt_height"]),
                            smoothing_length=as.numeric(params["smoothing"]),
                            useRandomDecoyModel=TRUE)
  
  # res[,corr:=as.numeric(params["corr"])]
  # res[,window:=as.numeric(params["window"])]
  # res[,rt_height:=as.numeric(params["rt_height"])]
  # res[,smoothing_length:=as.numeric(params["smoothing"])]
  res[]
}
