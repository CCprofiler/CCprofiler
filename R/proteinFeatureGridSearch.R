#' generateRandomPepTraces
#' @description generateRandomPepTraces.
#' @param traces traces object of type peptide
#' @param append Logical, whether to append the generated decoy traces
#' to the original target traces
#' @param traces traces bject of type peptide
#' @return traces object of type peptide with random peptides per protein
#' @export
#' @examples 
#' ## Load example data
#' exampleTraces <- examplePeptideTraces
#' 
#' ## Generate the decoy traces object
#' decoyTraces <- generateRandomPepTraces(exampleTraces)
#' # Inspect the result
#' summary(decoyTraces)
#' #--------------------
#' ## Directly append decoy traces to target taces
#' tracesWithDecoys <- generateRandomPepTraces(exampleTraces,
#'                                             append = T)
#' # Inspect the result
#' summary(tracesWithDecoys)
#' 

generateRandomPepTraces <- function(traces,
                                    append = FALSE,
                                    rm_orig_decoys = TRUE){
  if (rm_orig_decoys) {
    n_decoys <- length(grep("^DECOY", traces$trace_annotation$protein_id))
    if ( n_decoys > 0){
      idx_decoys <- grep("^DECOY_",traces$trace_annotation$protein_id)
      traces$traces <- traces$traces[-idx_decoys]
      traces$trace_annotation<- traces$trace_annotation[-idx_decoys]
      message(n_decoys, " original decoys removed")
    } else {
      message("no original decoys contained/removed")
    }
  }
  peptides <- traces$trace_annotation$id
  set.seed(123)
  random_peptides <- sample(peptides)
  traces_random <- traces
  traces_random$trace_annotation$id <- paste0("DECOY_", random_peptides)
  traces_random$trace_annotation$protein_id <- paste0("DECOY_",traces_random$trace_annotation$protein_id)
  traces_random$traces$id <- paste0("DECOY_", traces_random$traces$id)
  setorder(traces_random$trace_annotation, id)
  if(append){
    traces_res <- traces
    traces_res$traces <- rbind(traces$traces, traces_random$traces)
    traces_res$trace_annotation <- rbind(traces$trace_annotation, traces_random$trace_annotation)
    return(traces_res)
  }else{
    return(traces_random)
  }
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
                            smoothing_length=as.numeric(params["smoothing"]),
                            useRandomDecoyModel=TRUE
                            )
  res[,corr:=as.numeric(params["corr"])]
  res[,window:=as.numeric(params["window"])]
  res[,rt_height:=as.numeric(params["rt_height"])]
  res[,smoothing_length:=as.numeric(params["smoothing"])]
  res[]
}
