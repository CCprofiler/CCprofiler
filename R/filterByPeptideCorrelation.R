# Due to: http://stackoverflow.com/questions/24501245/data-table-throws-object-not-found-error
.datatable.aware=TRUE


filterByPeptideCorrelation <- function(Traces, cutoff= 0.4 , Keep1pep = FALSE) {
  
  df <- as.data.frame(Traces$traces)
  label <- as.data.frame(Traces$trace_annotation) 
  df <- merge(df, label, by = 'id' , all.x = TRUE)
  
  if (Keep1pep) {
    df[is.na(df$SibPepCorr),]$SibPepCorr <- 1
  }
  
  
  df <- subset(df, df$SibPepCorr >= cutoff)
  
  traces <- subset(df, select = c(2:(ncol(df)-3) , 1))
  trace_type <- 'peptide'
  trace_annotation <- subset(df, select = c(1,(ncol(df)-2):ncol(df)))
  fraction_annotation <- Traces$fraction_annotation
  
  traces <- as.data.table(traces)
  trace_annotation <- as.data.table(trace_annotation)
  
  result <- list("traces" = traces,
                 "trace_type" = trace_type,
                 "trace_annotation" = trace_annotation,
                 "fraction_annotation" = fraction_annotation)
  class(result) <- "traces"
  names(result) <- c("traces", "trace_type", "trace_annotation", "fraction_annotation")
  
  return(result) #return the filtered Traces data
  
}
