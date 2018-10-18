#' Summary of alternative peptide sequences
#' @description Select most abndant peptide for each group of peptides with same
#' genomic star position. These are e.g. peptides with missed cleavages.
#' @param traces An object of type traces, trace_type must be peptide.
#' @param topN Numeric integer, specifying the number of peptides to sum for protein quantification. Default is 2.
#' @return An object of type traces, trace_type is peptide.
#' @export
summarizeAlternativePeptideSequences <- function(traces,
                                  topN = 1,position="PeptidePositionStart"){
  UseMethod("summarizeAlternativePeptideSequences", traces)
}

#' @describeIn summarizeAlternativePeptideSequences
#' @export
summarizeAlternativePeptideSequences.traces <- function(traces,topN=1,position="PeptidePositionStart"){
  .tracesTest(traces, "peptide")
  if (! position %in% names(traces$trace_annotation)) {
    stop("This function is only available for traces annotated by the
    relative protein sequence position (annotateRelativePepPos).")
  }
  traces$trace_annotation[,protein_id_original := protein_id]
  traces$trace_annotation[,protein_id := paste0(protein_id,"_",get(position))]
  traces_summed <- proteinQuantification.traces(traces,
                                  topN = topN,
                                  keep_less = TRUE,
                                  rm_decoys = TRUE,
                                  use_sibPepCorr = FALSE,
                                  use_repPepCorr = FALSE,
                                  full_intersect_only = FALSE,
                                  quantLevel = "protein_id")
  traces_summed_new <- copy(traces_summed)
  ann <- subset(traces_summed_new$trace_annotation,select=c("protein_id","quant_peptides_used"))
  traces_summed_new$traces <- merge(traces_summed_new$traces,ann,by.x="id",by.y="protein_id")
  traces_summed_new$traces[,id:=NULL]
  setnames(traces_summed_new$traces,"quant_peptides_used","id")
  setorder(traces_summed_new$traces, -id)

  traces_summed_new$trace_annotation[,id:=quant_peptides_used]
  traces_summed_new$trace_annotation[,quant_peptides_used:=NULL]
  traces_summed_new$trace_annotation[,protein_id:=protein_id_original]
  traces_summed_new$trace_annotation[,protein_id_original:=NULL]
  setnames(traces_summed_new$trace_annotation,"n_peptides","n_cleaved_peptides")
  setorder(traces_summed_new$trace_annotation, -id)

  traces_summed_new$trace_type="peptide"

  .tracesTest(traces_summed_new, "peptide")

  return(traces_summed_new)
}

#' @describeIn summarizeAlternativePeptideSequences
#' @export
summarizeAlternativePeptideSequences.tracesList <- function(tracesList,topN=1){
  .tracesListTest(tracesList, type = "peptide",position="PeptidePositionStart")
  res_list <- lapply(tracesList, summarizeAlternativePeptideSequences.traces,
                topN=topN,position=position)
  class(res) <- "tracesList"
  .tracesListTest(res)
  return(res)
}
