#' Summary of alternative peptide sequences
#' @description Select most abndant peptide for each group of peptides with same
#' genomic star position. These are e.g. peptides with missed cleavages.
#' @param traces An object of type traces, trace_type must be peptide.
#' @param topN Numeric integer, specifying the number of peptides to sum for protein quantification. Default is 2.
#' @param verbose Logical if warning messages should be printed. Default is TRUE.
#' @return An object of type traces, trace_type is peptide.
#' @export
summarizeAlternativePeptideSequences <- function(traces,
                                  topN = 1,
                                  start="PeptidePositionStart",
                                  end="PeptidePositionEnd",
                                  verbose = TRUE){
  UseMethod("summarizeAlternativePeptideSequences", traces)
}

#' @describeIn summarizeAlternativePeptideSequences Select most abndant
#' peptide for each group of peptides with same
#' genomic star position. These are e.g. peptides with missed cleavages.
#' @export
summarizeAlternativePeptideSequences.traces <- function(traces,topN=1,
                                                        start="PeptidePositionStart",
                                                        end="PeptidePositionEnd",
                                                        verbose = TRUE){
  .tracesTest(traces, "peptide")
  if (! start %in% names(traces$trace_annotation)) {
    stop("This function is only available for traces annotated by the
    relative protein sequence position (annotateRelativePepPos).")
  }
  if (! end %in% names(traces$trace_annotation)) {
    stop("This function is only available for traces annotated by the
    relative protein sequence position (annotateRelativePepPos).")
  }
  annotation_old <- copy(traces$trace_annotation)
  traces_old <- copy(traces)
  traces_old$trace_annotation[,protein_id_original := protein_id]
  
  traces_old$trace_annotation[, new_id := "new"]
  
  proteins <- unique(traces_old$trace_annotation$protein_id)
  for(p in proteins){
    #print(p)
    p_ann <- subset(traces_old$trace_annotation, protein_id==p)
    peptides <- unique(p_ann$id)
    pep_seq <- list()
    for (pep in peptides) {
      pep_ann <- subset(p_ann, id==pep)
      coPeps <- subset(p_ann, ((p_ann[[start]] >= pep_ann[[start]]) & (p_ann[[start]] <= pep_ann[[end]])) | 
                         (p_ann[[end]] >= pep_ann[[start]]) & (p_ann[[end]] <= pep_ann[[end]])
      )$id
      pep_seq[[pep]] <- coPeps
    }
    pep_seq <- unique(pep_seq)
    #all_peps <- list()
    for (pep in peptides) {
      #all_peps[[pep]] <- paste0(as.character(sort(unique(unlist(lapply(pep_seq, function(x){if(pep %in% x)x}))))), collapse=";")
      traces_old$trace_annotation[(protein_id==p) & (id==pep)]$new_id <- paste0(as.character(sort(unique(unlist(lapply(pep_seq, function(x){if(pep %in% x)x}))))), collapse=";")
    }
  }
  
  traces_old$trace_annotation[,protein_id := paste0(protein_id,"_",new_id)]
  
  traces_summed <- proteinQuantification(traces_old,   #.traces
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
  
  setorder(traces_summed_new$trace_annotation, -id)
  
  traces_summed_new$trace_type="peptide"
  
  .tracesTest(traces_summed_new, "peptide")
  
  return(traces_summed_new)
}

#' @describeIn summarizeAlternativePeptideSequences Select most abndant
#' peptide for each group of peptides with same
#' @export
summarizeAlternativePeptideSequences.tracesList <- function(traces,topN=1,
                                                        start="PeptidePositionStart",
                                                        end="PeptidePositionEnd",
                                                        verbose = TRUE){
  .tracesListTest(tracesList, type = "peptide")
  res_list <- lapply(tracesList, summarizeAlternativePeptideSequences.traces,
                     topN=topN,
                     start=start,
                     end=end,
                     verbose = verbose)
  class(res) <- "tracesList"
  .tracesListTest(res)
  return(res)
}


#' Annotate peptides with sequence information
#' @description Annotate peptides with sequence information.
#' @param traces An object of type traces, trace_type must be peptide.
#' @param fasta_file path to a fasta file.
#' @return An object of type traces, trace_type is peptide.
#' @importFrom seqinr words.pos
#' @importFrom Biostrings readAAStringSet
#' @importFrom Biostrings toString
#' @export
annotatePeptideSequences <- function(traces,fasta_file){
  UseMethod("annotatePeptideSequences", traces)
}

#' @describeIn annotatePeptideSequences Annotate peptides with sequence information
#' @export
annotatePeptideSequences.traces <- function(traces,fasta_file){
  # Depends on:
  # library("seqinr")
  # library("Biostrings")
  fasta <- readAAStringSet(fasta_file)
  names(fasta) = gsub(".*\\|(.*?)\\|.*", "\\1", names(fasta))
  print("Following proteins are not in the provided fasta:")
  traces$trace_annotation[,PeptidePositionStart := getPepStartSite(id, protein_id, fasta), by=c("id","protein_id")]
  traces$trace_annotation[,PeptidePositionEnd := PeptidePositionStart+nchar(gsub("\\(UniMod:[0-9]+\\)","",id))-1, by=c("id","protein_id")]
  return(traces)
}
  
#' @describeIn annotatePeptideSequences Annotate peptides with sequence information
#' @export
annotatePeptideSequences.tracesList <- function(traces,fasta_file){
  .tracesListTest(traces, type = "peptide")
  res <- lapply(traces, annotatePeptideSequences.traces,
                     fasta_file=fasta_file)
  class(res) <- "tracesList"
  .tracesListTest(res)
  return(res)
}

# helper function
getPepStartSite <- function(pep, prot, fasta){
  if (prot %in% names(fasta)){
    return(words.pos(gsub("\\(UniMod:[0-9]+\\)","",pep), toString(fasta[prot]))[1])
  } else {
    print(prot)
    return(NA)
  }
}
