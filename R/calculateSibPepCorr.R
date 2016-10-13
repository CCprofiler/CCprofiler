#' Calculate sibling peptide correlation in traces.object of type peptide
#' @param traces.obj An object of type \code{traces.obj}.
#' @param plot logical TRUE or FALSE 
#' @param PDF logical TRUE or FALSE 
#' @param FFT numeric fraction of false targets (pi0) 
#' @return traces.obj An object of type \code{traces.obj}.
#' @export

calculateSibPepCorr <- function(pepTraces,
                                plot = TRUE,
                                PDF = FALSE,
                                FFT = 0.5) {
  
  if (pepTraces$trace_type != "peptide"){
    stop("Sibling peptide correlation can only be calculated on traces of type peptide")
  }
  
  quantdata <- getIntensityMatrix.traces(pepTraces)
  proteins <- unique(pepTraces$trace_annotation$protein_id)
  nproteins <- length(proteins)
  
  passed <- logical(length = nrow(pepTraces$trace_annotation))
  SibPepCorr <- numeric(length = nrow(pepTraces$trace_annotation))
  
  
  for (i in 1:nproteins){
    
    message(paste("PROCESSED", i, "of", nproteins, "proteins"))
    #indexpos <- proteins[i] == data$protein_id
    indexpos <- proteins[i] == pepTraces$trace_annotation$protein_id
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
  
  # Result <- cbind(labels, SibPepCorr, quantdata)
  
#   if (!Keep1pep){
#     Result <- subset(Result, !is.na(SibPepCorr))
#   }
#   
  #return the filtered Traces data
#   traces <- subset(Result, select = c(4:ncol(Result),1))
#   trace_type <- 'peptide'
  pepTraces$trace_annotation$SibPepCorr <- SibPepCorr
  # output plots
  if (plot){
    if (PDF){
      pdf("SibPepCorrFilter_plot.pdf")
    }
    plot.SibPepCorrDensities(pepTraces)
    ROC.SibPepCorr(pepTraces,FFT = FFT)
    if (PDF){
      dev.off()
    }
  }
  return(pepTraces)
}



