# Due to: http://stackoverflow.com/questions/24501245/data-table-throws-object-not-found-error
.datatable.aware=TRUE

#' calculateSibPepCorr
#' @description Calculate sibling peptide correlation in traces.object of type peptide.
#' @import data.table
#' @param traces An object of type traces.
#' @param plot logical,wether to print SibPepCorr density plot to R console.
#' @param PDF logical, wether to print SibPepCorr density plot to a PDF file.
#' @param name Character string with name of the plot, only used if \code{PDF=TRUE}.
#' PDF file is saved under name.pdf. Default is "SibPepCorr_densityplot".
#' @return An object of type traces with added SibPepCorr column.
#' @export
#' @example 
#' ## Load example data
#'  tracesRaw <- examplePeptideTracesUnannotated
#'  
#'  ## Calculate the SibPepCorr of every peptide
#'  tracesRawSpc <- calculateSibPepCorr(traces = tracesRaw,
#'                                       plot = TRUE)
#'                    
#'  head(tracesRawSpc$trace_annotation)                                     
#'  
#' 
calculateSibPepCorr <- function(traces,
                                plot = TRUE,
                                PDF = FALSE,
                                name = "SibPepCorr_densityplot")
  {

  # Check input type
  if (traces$trace_type != "peptide"){
    stop("Sibling peptide correlation can only be calculated on traces of type peptide")
  }

  # prepare data
  quantdata <- getIntensityMatrix(traces)
  proteins <- unique(traces$trace_annotation$protein_id)
  nproteins <- length(proteins)

  # calculation of SibPepCorr
  passed <- logical(length = nproteins)
  SibPepCorr <- numeric(length = nproteins)
  for (i in 1:nproteins){

    #message(paste("PROCESSED", i, "of", nproteins, "proteins"))
    #indexpos <- proteins[i] == data$protein_id
    indexpos <- proteins[i] == traces$trace_annotation$protein_id
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

  traces$trace_annotation$SibPepCorr <- SibPepCorr

  # output plot
  if (plot){
    if (PDF){
      name <- gsub("$|\\.pdf$", ".pdf", name)
      pdf("SibPepCorr_densityplot.pdf")
    }
    plotSibPepCorrDensities(traces)
    if (PDF){
      dev.off()
    }
  }
  return(traces)
}
