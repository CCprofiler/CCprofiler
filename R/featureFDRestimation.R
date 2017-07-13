#' Estimate decoy based featureFinding FDR.
#' @description Estimate FDR statistics based on a decoy model after feature finding.
#' @param features data.table containing filtered complex feature results.
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
#' ## estimate the FDR
#' estimateDecoyFDR(complexFeatures, FFT = 1)
#' 
#' #-----------------------------
#' ## Protein Features
#' #-----------------------------
#' 
#' ## Load example data
#' proteinFeatures <-exampleProteinFeatures
#' 
#' ## estimate the FDR
#' estimateDecoyFDR(proteinFeatures, FFT = 1)
#' 

estimateDecoyFDR <- function(features,
                             FFT = 1){
  
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
    message("No decoy features found. Please check if decoys were available in the complex hypotheses.")
    FDR <- 0
  } else {
    # estimate FDR based on detected decoys
    FDR <- (FFT *length(decoys)) / length(targets)
  }
  # estimate number of true positive complexes by correcting with the estimated FDR
  P <- length(targets)
  TP <- length(targets)*(1-FDR)
  D <- length(decoys)
  #output
  cols <- names(features)
  if("corr" %in% cols & "window" %in% cols & "rt_height" %in% cols & "smoothing_length" %in% cols& "peak_corr_cutoff" %in% cols& "completeness_cutoff" %in% cols& "n_subunits_cutoff" %in% cols){
    list(FDR = FDR, TP = TP, P = P, D = D,
         corr = unique(features$corr),
         window = unique(features$window),
         rt_height = unique(features$rt_height),
         smoothing_length = unique(features$smoothing_length),
         peak_corr_cutoff = unique(features$peak_corr_cutoff),
         completeness_cutoff = unique(features$completeness_cutoff),
         n_subunits_cutoff = unique(features$n_subunits_cutoff))
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
