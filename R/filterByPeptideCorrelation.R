# Due to: http://stackoverflow.com/questions/24501245/data-table-throws-object-not-found-error
.datatable.aware=TRUE


filterByPeptideCorrelation <- function(Traces, cutoff= 0.4 , Keep1pep = FALSE) {
  
  df <- as.data.frame(Traces$traces)
  label <- as.data.frame(Traces$traces_annotation) 
  df <- merge(df, label, by = 'id' , all.x = TRUE)
  
  if (Keep1pep) {
    df[is.na(df$SibPepCorr),]$SibPepCorr <- 1
  }
  
  
  df <- subset(df, df$SibPepCorr >= cutoff)
  
  traces <- subset(df, select = c(2:(ncol(df)-3) , 1))
  traces_type <- 'peptide'
  traces_annotation <- subset(df, select = c(1,(ncol(df)-2):ncol(df)))
  fraction_annotation <- Traces$fraction_annotation
  
  traces <- as.data.table(traces)
  traces_annotation <- as.data.table(traces_annotation)
  
  result <- list("traces" = traces,
                 "traces_type" = traces_type,
                 "traces_annotation" = traces_annotation,
                 "fraction_annotation" = fraction_annotation)
  class(result) <- "Traces"
  names(result) <- c("traces", "traces_type", "traces_annotation", "fraction_annotation")
  
  return(result) #return the filtered Traces data
  
}
