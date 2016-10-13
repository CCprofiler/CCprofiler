# Due to: http://stackoverflow.com/questions/24501245/data-table-throws-object-not-found-error
.datatable.aware=TRUE


filterByPeptideCorrelation <- function(Traces, cutoff= 0.4 , Keep1pep = FALSE) {
  
  if (Keep1pep) {
    Traces$trace_annotation[is.na(SibPepCorr),SibPepCorr:=1]
  }
  
  peptides_to_keep <- Traces$trace_annotation$id[which(Traces$trace_annotation$SibPepCorr >= cutoff)]
  Traces$traces <- subset(Traces$traces,id %in% peptides_to_keep)
  Traces$trace_annotation<- subset(Traces$trace_annotation,id %in% peptides_to_keep)

  return(Traces) #return the filtered Traces data
  
}
