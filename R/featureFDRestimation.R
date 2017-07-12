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
