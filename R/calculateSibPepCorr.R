#' Calculate sibling peptide correlation in traces.object of type peptide
#' @param traces.obj An object of type \code{traces.obj}.
#' @param plot logical TRUE or FALSE 
#' @param PDF logical TRUE or FALSE 
#' @param FFT numeric fraction of false targets (pi0) 
#' @return traces.obj An object of type \code{traces.obj}.
#' @export

calculateSibPepCorr <- function(Traces,
                                plot = TRUE,
                                PDF = FALSE)
  {
  
  # Check input type
  if (Traces$trace_type != "peptide"){
    stop("Sibling peptide correlation can only be calculated on traces of type peptide")
  }
  
  # prepare data
  quantdata <- getIntensityMatrix.traces(Traces)
  proteins <- unique(Traces$trace_annotation$protein_id)
  nproteins <- length(proteins)
  
  # calculation of SibPepCorr
  passed <- logical(length = nproteins)
  SibPepCorr <- numeric(length = nproteins)
  for (i in 1:nproteins){
    
    message(paste("PROCESSED", i, "of", nproteins, "proteins"))
    #indexpos <- proteins[i] == data$protein_id
    indexpos <- proteins[i] == Traces$trace_annotation$protein_id
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
    } else { #if only one peptide, check with the option
      SibPepCorr[indexpos] <- NA
    }
  }
  
  Traces$trace_annotation$SibPepCorr <- SibPepCorr
  
  # output plot
  if (plot){
    if (PDF){
      pdf("SibPepCorr_densityplot.pdf")
    }
    plot.SibPepCorrDensities(Traces)
    if (PDF){
      dev.off()
    }
  }
  return(Traces)
}



