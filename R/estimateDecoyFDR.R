#' Estimate statistics based on decoy model.
#' @description Estimate statistics based on decoy model.
#' @param complex_features data.table containing filtered complex feature results.
#' @return List with stats
#' @export
estimateDecoyFDR <- function(complex_features){
  complexes <- complex_features$complex_id
  decoys = complexes[grep("DECOY",complexes)]
  targets = complexes[!complexes %in% decoys]
  if(!length(decoys)>0){
    message("No decoy features found. Please check if decoys were available in the complex hypotheses.")
    FDR = 0
  } else {
    # estimate FDR based on detected decoys
    FDR = length(decoys)/length(targets)
  }
  # estimate number of true positive complexes by correcting with the estimated FDR
  P = length(targets)
  TP = length(targets)*(1-FDR)
  list(FDR=FDR,TP=TP,P=P)
}
