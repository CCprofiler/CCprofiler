# Due to: http://stackoverflow.com/questions/24501245/data-table-throws-object-not-found-error
# .datatable.aware=TRUE

#' Import peptide profiles from an OpenSWATH experiment.
#' @name filterByPeptideCorrelation
#' @import data.table
#' @description filters peptide elution profiles based on sibling peptide correlation. 
#' @param Traces Object of class Traces containing the quantitative MS data table.
#' @param min_pairwise_corr Minimum pairwise pearson correlation necessary to retain the data for a protein with 2 peptides.
#' @param average_corr_cutoff For proteins with >= 3 peptides, retain all in the data whose average correlation with all siblings is
#'       how good, defaults to ">= median" but can be replaces with a numeric cutoff (e.g. 0.6).
#'        
#' @return An object of class Traces (peptide level) containing table.wide, table.long and id} that can be processed with the herein contained functions.
#'         
#' @export

filterByPeptideCorrelation <- function(Traces, min_pairwise_corr = 0.8, average_corr_cutoff = ">= median"){
  
  data <- Traces$traces.wide
  quantdata <- as.matrix(data[,3:ncol(data), with = FALSE])
  labels <- data[,1:2]
  rownames(quantdata) <- data$peptide_id
  proteins <- unique(data$protein_id)
  nproteins <- length(proteins)
  passed <- logical(length = nrow(data))
  for (i in 1:nproteins){
    message(paste("PROCESSED", i, "of", nproteins, "proteins"))
    indexpos <- proteins[i] == data$protein_id
    df <- quantdata[indexpos,]
    df_cor <- cor(t(df))
    if (isTRUE(nrow(df) == 2)){ # If there's only 2 peptides, request minpaircorr
      paircorr <- mean(rowMeans(df_cor))
      if (isTRUE(paircorr >= min_pairwise_corr)) {
        passed[indexpos] <- TRUE
      } else{
        passed[indexpos] <- FALSE
      }
    } else if (isTRUE(nrow(df) >= 3)) { 
      # If there's 3 or more peps, keep those correlated
      # equal or better than median
      avgcorr <- rowMeans(df_cor)
      if (average_corr_cutoff == ">= median"){
        survivors <- names(avgcorr[avgcorr >= median(avgcorr) ] )
      }
      if (class(average_corr_cutoff) == "numeric"){
        survivors <- names(avgcorr[avgcorr >= average_corr_cutoff ] )
      }
      indexpos <- data$peptide_id %in% survivors
      passed[indexpos] <- TRUE
    } else { # if only one peptide, do not keep
      passed[indexpos] <- FALSE
    }
  }
  Result <- list(traces.wide = data[passed,], ids = Traces$ids[passed,]) 
  class(Result) <- "Traces"
  Result #return the filtered Traces data
}
