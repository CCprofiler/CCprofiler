

# For some pretty weird reason
.datatable.aware=TRUE


CalculateSibPepCorr <- function(pepTraces, Keep1pep = FALSE , Density_plot = TRUE , ROC_curve = TRUE , FFT = 0.5) {
  
  labels <- as.data.frame(pepTraces$traces_annotation)
  labels <- subset(labels, select = c('FullPeptideName','ProteinName'))
  data <- as.data.frame(pepTraces$traces)
  data <- merge(data, labels, by.x = 'id', by.y = 'FullPeptideName',all.x = TRUE)
  data <- as.data.frame(data[,c((2:(ncol(data)-1)),1,ncol(data))])
  labels <- subset (data, select = c("id", "ProteinName"))
  quantdata <- as.matrix(data[,1:(ncol(data)-2)])
  rownames(quantdata) <- data$id
  proteins <- unique(labels$ProteinName)
  nproteins <- length(proteins)
  passed <- logical(length = nrow(data))
  SibPepCorr <- numeric(length = nrow(data))
  
  
  for (i in 1:nproteins){
    
    message(paste("PROCESSED", i, "of", nproteins, "proteins"))
    indexpos <- proteins[i] == data$ProteinName
    df <- quantdata[indexpos,]
    class(df) <- 'numeric'
    df_cor <- cor(t(df))
    npep = nrow(df)
    
    if (isTRUE(npep == 2)){ # If there's only 2 peptides, request minpaircorr
      paircorr <- df_cor[1,2]
      SibPepCorr[indexpos] <- paircorr
    } else if (isTRUE(npep >= 3)) { 
      # If there's 3 or more peps
      sibcorrs <- sapply(1:nrow(df_cor), function(x){(sum(df_cor[x,])-1)/(nrow(df_cor)-1)})
      SibPepCorr[indexpos] <- sibcorrs
    } else {       #if only one peptide, check with the option
      SibPepCorr[indexpos] <- NA
    }
  }
  
  
  Result <- cbind(labels, SibPepCorr, quantdata)
  
  if (!Keep1pep){
    Result <- subset(Result, !is.na(SibPepCorr))
  }
  

  #return the filtered Traces data
  traces <- subset(Result, select = c(4:ncol(Result),1))
  traces_type <- 'peptide'
  traces_annotation <- subset(Result, select = c(1:3))
  traces_annotation$FullPeptideName <- traces_annotation$id
  fraction_annotation <- pepTraces$fraction_annotation
  
  results <- list("traces" = traces,
                 "traces_type" = traces_type,
                 "traces_annotation" = traces_annotation,
                 "fraction_annotation" = fraction_annotation)
  class(results) <- "Traces"
  names(results) <- c("traces", "traces_type", "traces_annotation", "fraction_annotation")
  
  if (Density_plot == TRUE) {
    plot.SibPepCorrDensities(results)
  }
  
  if (ROC_curve == TRUE) {
    ROC.SibPepCorr(results,FFT = FFT)
  }
  
  traces <- as.data.table(traces)
  traces_annotation <- as.data.table(traces_annotation)
  
  results <- list("traces" = traces,
                  "traces_type" = traces_type,
                  "traces_annotation" = traces_annotation,
                  "fraction_annotation" = fraction_annotation)
  class(results) <- "Traces"
  names(results) <- c("traces", "traces_type", "traces_annotation", "fraction_annotation")
  
  return(results)
}



