# Due to: http://stackoverflow.com/questions/24501245/data-table-throws-object-not-found-error
.datatable.aware=TRUE

#' Protein Quantification
#' @description Calculate protein quantities basen on the topN peptide intensities.
#' @import data.table
#' @param traces An object of type traces, trace_type must be peptide.
#' @param topN Numeric integer, specifying the number of peptides to sum for protein quantification. Default is 2.
#' @param keep_less Logical, specifying whether proteins with less than topN peptides
#'  should be kept in the data (This may result in some protein intensities being calculated
#'  as the sum of fewer peptide intensities than others. Only use withh caution.), default is \code{FALSE}.
#' @param rm_decoys Logical, specifying whether decoys should be kept.
#'  The decoys have only limited use on the protein level, default is \code{TRUE}.
#' @return An object of type traces, trace_type is protein.
#' @export
#' @example
#' ## Load example data
#' pepTraces <- examplePeptideTracesFiltered
#' 
#' ## Sum the intensities of the top 2 peptides to get protein intensities
#' protTraces <- proteinQuantification(pepTraces,
#'                                      topN = 2)
#' ## Check the result
#' summary(pepTraces)
#' summary(protTraces)
#' 
#' # ProteinTraces annotation
#' head(protTraces$trace_annotation)
#' # The protein_id column from the peptide traces object becomes the new id column of the protein traces object.
#' # The last 2 columns indicate how many peptides could be observed and which were summed for quantification.
#' 
#' 
proteinQuantification <- function(traces,
                                  topN = 2,
                                  keep_less = FALSE,
                                  rm_decoys = TRUE){
  
  ## Check if it's a peptide level table
  .tracesTest(traces, "peptide")
  
  ## remove decoys
  if (rm_decoys) {
    n_decoys <- length(grep("^DECOY", traces$trace_annotation$protein_id))
    if ( n_decoys > 0){
      idx_decoys <- grep("^DECOY_",traces$trace_annotation$protein_id)
      traces$traces <- traces$traces[-idx_decoys]
      traces$trace_annotation<- traces$trace_annotation[-idx_decoys]
      message(n_decoys, " decoys removed")
    } else {
      message("no decoys contained/removed")
    }
  }
  
  ## Extract wide table for calculations
  peptideTracesTable <- data.table(protein_id = traces$trace_annotation$protein_id,
                                   peptide_id = traces$trace_annotation$id,
                                   subset(traces$traces, select =-id))
  
  # Calculations in long format - sum the topN peptides per protein
  peptideTracesLong <- melt(peptideTracesTable,
                            id.vars = c("protein_id", "peptide_id"),
                            variable.name = "fraction_number",
                            value.name = "intensity")
  peptideTracesLong[, intensity:=as.numeric(intensity)]
  peptideTracesLong[, peptide_intensity:=sum(intensity), peptide_id]
  peptideTracesLong[, n_peptides:=length(unique(peptide_id)), protein_id]
  ## the ties.method makes sure how to deal with peptides of identical intensity: "first" keeps the order of occurence
  peptideTracesLong[, peptide_rank:=rank(-peptide_intensity[1:n_peptides[1]],ties.method = "first"), protein_id]
  peptideTracesLong <- peptideTracesLong[peptide_rank <= topN]
  ## collect information which peptides were used for quantification
  peptideTracesLong[, quant_peptides_used:=paste(unique(peptide_id), collapse = ","), protein_id]
  
  ## if not wanted, kick out those with less than N peptides
  if(!keep_less){
    peptideTracesLong <- peptideTracesLong[n_peptides >= topN]
  }
  
  ## Sum peptides to protein level (wide) traces table
  peptideTracesTopNsumWide <- as.data.table(cast(peptideTracesLong,
                                                 protein_id ~ fraction_number,
                                                 value = "intensity",
                                                 fun.aggregate = sum))
  
  ## move id column to end to ensure correct quant value index
  peptideTracesTopNsumWide[, id:=protein_id]
  peptideTracesTopNsumWide <- subset(peptideTracesTopNsumWide, select=-protein_id)
  setorder(peptideTracesTopNsumWide, -id)
  
  ## assemble updated, protein-level trace_annotation table
  oldAnnotationPeptidelevel = subset(traces$trace_annotation, select =-id)
  
  
  if ("SibPepCorr" %in% names(oldAnnotationPeptidelevel)){
    oldAnnotationPeptidelevel[, SibPepCorr_protein_mean:=mean(SibPepCorr), protein_id]
    ## removal of peptide level SibPepCorr necessary after merge.data.table update 2016-11
    oldAnnotationPeptidelevel = unique(subset(oldAnnotationPeptidelevel, select =-SibPepCorr))
  }
  new_annotation = unique(peptideTracesLong[,.(protein_id, n_peptides, quant_peptides_used)])
  peptideTracesTopNsumWideAnnotation <- merge(oldAnnotationPeptidelevel,
                                              new_annotation,
                                              by = "protein_id", all.x = FALSE, all.y = TRUE)
  
  peptideTracesTopNsumWideAnnotation <- unique(peptideTracesTopNsumWideAnnotation)
  peptideTracesTopNsumWideAnnotation[, id:=protein_id]
  ## Order traces alphabetically according to their id and 
  setorder(peptideTracesTopNsumWideAnnotation, -id)
  setcolorder(peptideTracesTopNsumWideAnnotation, c("id", names(peptideTracesTopNsumWideAnnotation[,!"id", with = F])))
  ## assemble result protein level traces object
  traces$traces <- peptideTracesTopNsumWide
  traces$trace_annotation <- peptideTracesTopNsumWideAnnotation
  traces$trace_type <- "protein"
  
  ## Test and return
  .tracesTest(traces, type = "protein")
  return(traces)
}
