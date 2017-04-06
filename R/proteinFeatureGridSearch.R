#' generateRandomPepTraces
#' @description generateRandomPepTraces.
#' @param traces traces bject of type peptide
#' @return traces object of type peptide with random peptides per protein
#' @export
generateRandomPepTraces <- function(traces){
  peptides = traces$trace_annotation$id
  set.seed(123)
  random_peptides = sample(peptides)
  traces_random = traces
  traces_random$trace_annotation$id=random_peptides
  traces_random$trace_annotation$protein_id <- paste0("DECOY_",traces_random$trace_annotation$protein_id)
  traces_random
}

#' getBestProteinFeatures
#' @description getBestProteinFeatures.
#' @param res data.table with protein feature finding results
#' @return List with only one best feature per protein
#' @export
getBestProteinFeatures <- function(res){
  res_best <- unique(res,by="protein_id")
  res_best
}

#' Perform Protein feature grid search
#' @description Perform Protein feature grid search.
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
performProteinGridSearch <- function(traces,
                                    calibration,
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
  # setting a seed is absolutely crutial to ensure reproducible results!!!!!!!!!!!!!!!!!!!
  clusterSetRNGStream(cl,123)
  doSNOW::registerDoSNOW(cl)
  clusterEvalQ(cl,library(SECprofiler,data.table))
  clusterExport(cl, list("generateRandomPepTraces"))
  data <- list()
  data <- c(data,parRapply(cl,parameter_grid,FUN=runGridProteinFeatureFinding,traces=traces,calibration=calibration))
  stopCluster(cl)
  data
}

runGridProteinFeatureFinding <- function(params,traces,calibration) {
  res = findProteinFeatures(traces=traces,
                            calibration = calibration,
                            corr_cutoff = as.numeric(params["corr"]),
                            window_size = as.numeric(params["window"]),
                            parallelized = FALSE,
                            collapse_method="apex_only",
                            perturb_cutoff="1%",
                            rt_height=as.numeric(params["rt_height"]),
                            smoothing_length=as.numeric(params["smoothing"])
                            )
  #generate random peptide traces
  traces_random <- generateRandomPepTraces(traces)
  res_random = findProteinFeatures(traces=traces_random,
                            calibration = calibration,
                            corr_cutoff = as.numeric(params["corr"]),
                            window_size = as.numeric(params["window"]),
                            parallelized = FALSE,
                            collapse_method="apex_only",
                            perturb_cutoff="1%",
                            rt_height=as.numeric(params["rt_height"]),
                            smoothing_length=as.numeric(params["smoothing"])
                            )
  res_all <- rbind(res,res_random)
  res_all[,corr:=as.numeric(params["corr"])]
  res_all[,window:=as.numeric(params["window"])]
  res_all[,rt_height:=as.numeric(params["rt_height"])]
  res_all[,smoothing_length:=as.numeric(params["smoothing"])]
  res_all[]
}
