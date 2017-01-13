proteinQuantification <- function(Traces, topN = 2, keep_less = FALSE, remove.decoys = TRUE){
  # Check if it's a peptide level table
  if(Traces$trace_type != "peptide"){
    stop("The input object is not a peptide level traces object but", Traces$trace_type)
  }

  #remove decoys
  if (remove.decoys) {
    n_decoys <- length(grep("^DECOY", Traces$trace_annotation$protein_id))
    if ( n_decoys > 0){
      idx_decoys <- grep("^DECOY_",Traces$trace_annotation$protein_id)
      Traces$traces <- Traces$traces[-idx_decoys]
      Traces$trace_annotation<- Traces$trace_annotation[-idx_decoys]
      message(n_decoys, " decoys removed")
    } else {
      message("no decoys contained/removed")
    }
  }

  # Extract wide table for calculations
  peptideTracesTable <- data.table(protein_id = Traces$trace_annotation$protein_id,
                                   peptide_id = Traces$trace_annotation$id,
                                   subset(Traces$traces, select =-id))

  # Calculations in long format - sum the topN peptides per protein
  peptideTraces.long <- melt(peptideTracesTable,
                             id.vars = c("protein_id", "peptide_id"),
                             variable.name = "fraction_number",
                             value.name = "intensity")
  peptideTraces.long[, intensity:=as.numeric(intensity)]
  peptideTraces.long[, peptide_intensity:=sum(intensity), peptide_id]
  peptideTraces.long[, n_peptides:=length(unique(peptide_id)), protein_id]
  # the ties.method makes sure how to deal with peptides of identical intensity: "first" keeps the order of occurence
  peptideTraces.long[, peptide_rank:=rank(-peptide_intensity[1:n_peptides[1]],ties.method = "first"), protein_id]
  peptideTraces.long <- peptideTraces.long[peptide_rank <= topN]
  # collect information which peptides were used for quantification
  peptideTraces.long[, quant_peptides_used:=paste(unique(peptide_id), collapse = ","), protein_id]

  # if not wanted, kick out those with less than N peptides
  if(!keep_less){
    peptideTraces.long <- peptideTraces.long[n_peptides >= topN]
  }

  # Sum peptides to protein level (wide) traces table
  peptideTraces.topNsum.wide <- as.data.table(cast(peptideTraces.long,
                                                   protein_id ~ fraction_number,
                                                   value = "intensity",
                                                   fun.aggregate = sum))

  # move id column to end to ensure correct quant value index
  peptideTraces.topNsum.wide[, id:=protein_id]
  peptideTraces.topNsum.wide <- subset(peptideTraces.topNsum.wide, select=-protein_id)
  setorder(peptideTraces.topNsum.wide, -id)

  ## assemble updated, protein-level trace_annotation table
  old_annotation_peptidelevel = subset(Traces$trace_annotation, select =-id)
  # removal of peptide level SibPepCorr necessary after merge.data.table update 2016-11
  if ("SibPepCorr" %in% names(old_annotation_peptidelevel)){
    old_annotation_peptidelevel[, SibPepCorr_protein_mean:=mean(SibPepCorr), protein_id]
	old_annotation_peptidelevel = unique(subset(old_annotation_peptidelevel, select =-SibPepCorr))
	}
  new_annotation = unique(peptideTraces.long[,.(protein_id, n_peptides, quant_peptides_used)])
  peptideTraces.topNsum.wide.annotation <- merge(old_annotation_peptidelevel,
                                                 new_annotation,
                                                 by = "protein_id", all.x = FALSE, all.y = TRUE)

  peptideTraces.topNsum.wide.annotation <- unique(peptideTraces.topNsum.wide.annotation)
  peptideTraces.topNsum.wide.annotation[, id:=protein_id]
  setorder(peptideTraces.topNsum.wide.annotation, -id)

  if (!all(peptideTraces.topNsum.wide.annotation$id == peptideTraces.topNsum.wide$id)){
    stop("traces$id != trace_annotation$id")
  }

  # assemble result protein level traces object
  Traces$traces <- peptideTraces.topNsum.wide
  Traces$trace_annotation <- peptideTraces.topNsum.wide.annotation
  Traces$trace_type <- "protein"
  return(Traces)
}
