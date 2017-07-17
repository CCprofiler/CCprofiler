#' Estimate decoy based featureFinding FDR.
#' @description Estimate FDR statistics based on a decoy model after feature finding.
#' @param features data.table containing filtered complex feature results.
#' @param FFT Numeric, fraction of false targets. Only useful in special cases, should not
#' be altered in a generic workflow. Defaults to \code{1}.
#' @param grid_search_params Character vector of columns to report with the statistics for the dataset.
#' Usually only useful for a grid search. Otherwhise the default should not be altered.
#' Default = \code{NULL}.
#' @param verbose Logical, wether to print messages. Default = \code{TRUE}.
#' @return List with stats
#' @export
#' @examples 
#' 
#' #-----------------------------
#' ## Complex Features
#' #-----------------------------
#' 
#' ## Load example data
#' complexFeatures <- exampleComplexFeatures
#' ## Filter feature table 
#' complexFeaturesFiltered <- filterFeatures(complexFeatures,
#'                                           min_peak_corr= 0.5,
#'                                           min_feature_completeness= 0.5) 
#'
#' ## estimate the FDR
#' estimateDecoyFDR(complexFeaturesFiltered)
#' 
#' #-----------------------------
#' ## Protein Features
#' #-----------------------------
#' 
#' ## Load example data
#' proteinFeatures <-exampleProteinFeatures
#' ## Filter feature table 
#' proteinFeaturesFiltered <- filterFeatures(proteinFeatures,
#'                                           min_peak_corr= 0.5,
#'                                           min_feature_completeness= 0.5) 
#' 
#' ## estimate the FDR
#' estimateDecoyFDR(proteinFeaturesFiltered)
#' 

estimateDecoyFDR <- function(features,
                             FFT = 1,
                             grid_search_params = NULL,
                             verbose = TRUE){
  
  if("complex_id" %in% names(features)){
    complexes <- features$complex_id
  } else if ("protein_id" %in% names(features)) {
    complexes <- features$protein_id
  } else {
    message("not a valid search result")
  }
  decoy_ind <- grepl("DECOY",complexes)
  decoys <- complexes[decoy_ind]
  targets <- complexes[!decoy_ind]
  if(!length(decoys)>0){
    if(verbose){
      message("No decoy features found. Please check if decoys were available in the complex hypotheses.")
    }
    FDR <- 0
  } else {
    # estimate FDR based on detected decoys
    FDR <- (FFT *length(decoys)) / length(targets)
  }
  # estimate number of true positive complexes by correcting with the estimated FDR
  P <- length(targets)
  TP <- length(targets)*(1-FDR)
  D <- length(decoys)
  
  # output
  if(!is.null(grid_search_params)){
    fdr_list <- list(FDR = FDR, TP = TP, P = P, D = D)
    param_list <- lapply(grid_search_params, function(x) as.numeric(unique(features[, x, with = F])))
    names(param_list) <- grid_search_params
    c(fdr_list, param_list)
  } else {
    list(FDR = FDR, TP = TP, P = P, D = D)
  }
}

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
#'                    
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
