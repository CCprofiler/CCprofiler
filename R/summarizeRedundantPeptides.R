#' Summary of alternative peptide sequences
#' @description Select most abndant peptide for each group of peptides with same
#' genomic star position. These are e.g. peptides with missed cleavages.
#' @param traces An object of type traces, trace_type must be peptide.
#' @param topN Numeric integer, specifying the number of peptides to sum for protein quantification. Default is 2.
#' @param verbose Logical if warning messages should be printed. Default is TRUE.
#' @return An object of type traces, trace_type is peptide.
#' @export
summarizeAlternativePeptideSequences <- function(traces,
                                  topN = 1,position="PeptidePositionStart",
                                  verbose = TRUE){
  UseMethod("summarizeAlternativePeptideSequences", traces)
}

#' @describeIn summarizeAlternativePeptideSequences Select most abndant
#' peptide for each group of peptides with same
#' genomic star position. These are e.g. peptides with missed cleavages.
#' @export
summarizeAlternativePeptideSequences.traces <- function(traces,topN=1,
  position="PeptidePositionStart",
  verbose = TRUE){
  .tracesTest(traces, "peptide")
  if (! position %in% names(traces$trace_annotation)) {
    stop("This function is only available for traces annotated by the
    relative protein sequence position (annotateRelativePepPos).")
  }
  annotation_old <- copy(traces$trace_annotation)
  traces_old <- copy(traces)
  traces_old$trace_annotation[,protein_id_original := protein_id]
  traces_old$trace_annotation[,protein_id := paste0(protein_id,"_",get(position))]
  traces_summed <- proteinQuantification.traces(traces_old,
                                  topN = topN,
                                  keep_less = TRUE,
                                  rm_decoys = TRUE,
                                  use_sibPepCorr = FALSE,
                                  use_repPepCorr = FALSE,
                                  full_intersect_only = FALSE,
                                  quantLevel = "protein_id",
                                  verbose = verbose)
  traces_summed_new <- copy(traces_summed)
  ann <- subset(traces_summed_new$trace_annotation,select=c("protein_id","quant_peptides_used"))
  traces_summed_new$traces <- merge(traces_summed_new$traces,ann,by.x="id",by.y="protein_id")
  traces_summed_new$traces[,id:=NULL]
  setnames(traces_summed_new$traces,"quant_peptides_used","id")
  setorder(traces_summed_new$traces, -id)

  traces_summed_new$trace_annotation <- subset(annotation_old,id %in% traces_summed_new$traces$id)
  #traces_summed_new$trace_annotation <- copy(traces$trace_annotation)
  #traces_summed_new$trace_annotation[,id:=quant_peptides_used]
  #traces_summed_new$trace_annotation[,quant_peptides_used:=NULL]
  #traces_summed_new$trace_annotation[,protein_id:=protein_id_original]
  #traces_summed_new$trace_annotation[,protein_id_original:=NULL]
  ##setnames(traces_summed_new$trace_annotation,"n_peptides","n_cleaved_peptides")
  #traces_summed_new$trace_annotation[,n_peptides:=NULL]
  setorder(traces_summed_new$trace_annotation, -id)

  traces_summed_new$trace_type="peptide"

  .tracesTest(traces_summed_new, "peptide")

  return(traces_summed_new)
}

#' @describeIn summarizeAlternativePeptideSequences Select most abndant
#' peptide for each group of peptides with same
#' @export
summarizeAlternativePeptideSequences.tracesList <- function(tracesList,topN=1,
  position="PeptidePositionStart",
  verbose = TRUE){
  .tracesListTest(tracesList, type = "peptide")
  res_list <- lapply(tracesList, summarizeAlternativePeptideSequences.traces,
                topN=topN,position=position,
                verbose = verbose)
  class(res) <- "tracesList"
  .tracesListTest(res)
  return(res)
}
