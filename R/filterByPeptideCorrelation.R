<<<<<<< HEAD
# Due to: http://stackoverflow.com/questions/24501245/data-table-throws-object-not-found-error
.datatable.aware=TRUE

#' Import peptide profiles from an OpenSWATH experiment.
#' @name filterByPeptideCorrelation
#' @import data.table
#' @description filters peptide elution profiles based on sibling peptide
#'     correlation. 
#' @param traces Object of class traces containing the quantitative MS data
#'     table.
#' @param min.pairwise.corr Minimum pairwise pearson correlation necessary
#'     to retain the data for a protein with 2 peptides.
#' @param average.corr.cutoff For proteins with >= 3 peptides, retain all
#'     in the data whose average correlation with all siblings is how good,
#'     defaults to '>= median' but can be replaces with a numeric cutoff
#'     (e.g. 0.6).
#' @return An object of class traces (peptide level) containing table.wide,
#'     table.long and id that can be processed with the herein contained
#'     functions.
#' @export
filterByPeptideCorrelation <- function(Traces,
                                       min.pairwise.corr=0.8,
                                       average.corr.cutoff='>= median') {
  quantdata <- getIntensityMatrix(Traces)
  labels <- merge(Traces$traces$id, Traces$traces_annotation, by = "id")
  proteins <- unique(labels$protein_id)
  nproteins <- length(proteins)
  passed <- logical(length=nrow(data))
  for (i in 1:nproteins) {
    message(paste('PROCESSED', i, 'of', nproteins, 'proteins'))
    indexpos <- proteins[i] == labels$protein_id
    df <- quantdata[indexpos, ]
    df_cor <- cor(t(df))
    if (isTRUE(nrow(df) == 2)){ # If there's only 2 peptides, request minpaircorr
      paircorr <- mean(df_cor[1,2])
      if (isTRUE(paircorr >= min.pairwise.corr)) {
        passed[indexpos] <- TRUE
      } else{
        passed[indexpos] <- FALSE
      }
    } else if (isTRUE(nrow(df) >= 3)) {
      # If there's 3 or more peps, keep those correlated
      # equal or better than cutoff
      avgcorr <- rowMeans(df_cor) - 1/nrow(df)
      if (average.corr.cutoff == '>= median') {
        survivors <- names(avgcorr[avgcorr >= median(avgcorr)])
      }
      if (class(average.corr.cutoff) == 'numeric') {
        survivors <- names(avgcorr[avgcorr >= average.corr.cutoff])
      }
      indexpos <- data$peptide_id %in% survivors
      passed[indexpos] <- TRUE
    } else { # if only one peptide, do not keep
      passed[indexpos] <- FALSE
    }
  }
  traces_filtered <- Traces$traces[passed,]
  traces_annotation_filtered <- Traces$traces_annotation[passed,]
  
  result <- Traces
  result$traces <- traces_filtered
  result$traces_annotation <- traces_annotation_filtered
  return(result) #return the filtered Traces data
}
=======
# Due to: http://stackoverflow.com/questions/24501245/data-table-throws-object-not-found-error
# .datatable.aware=TRUE

#' Import peptide profiles from an OpenSWATH experiment.
#' @name filterByPeptideCorrelation
#' @import data.table
#' @description filters peptide elution profiles based on sibling peptide
#'     correlation. 
#' @param traces Object of class traces containing the quantitative MS data
#'     table.
#' @param min.pairwise.corr Minimum pairwise pearson correlation necessary
#'     to retain the data for a protein with 2 peptides.
#' @param average.corr.cutoff For proteins with >= 3 peptides, retain all
#'     in the data whose average correlation with all siblings is how good,
#'     defaults to '>= median' but can be replaces with a numeric cutoff
#'     (e.g. 0.6).
#' @return An object of class traces (peptide level) containing table.wide,
#'     table.long and id that can be processed with the herein contained
#'     functions.
#' @export
filterByPeptideCorrelation <- function(traces,
                                       min.pairwise.corr=0.8,
                                       average.corr.cutoff='>= median') {
  data <- traces$traces.wide
  quantdata <- as.matrix(data[, 3:ncol(data), with=FALSE])
  labels <- data[, 1:2]
  rownames(quantdata) <- data$peptide_id
  proteins <- unique(data$protein_id)
  nproteins <- length(proteins)
  passed <- logical(length=nrow(data))
  for (i in 1:nproteins) {
    message(paste('PROCESSED', i, 'of', nproteins, 'proteins'))
    indexpos <- proteins[i] == data$protein_id
    df <- quantdata[indexpos, ]
    df_cor <- cor(t(df))
    if (isTRUE(nrow(df) == 2)){ # If there's only 2 peptides, request minpaircorr
      paircorr <- mean(rowMeans(df_cor))
      if (isTRUE(paircorr >= min.pairwise.corr)) {
        passed[indexpos] <- TRUE
      } else{
        passed[indexpos] <- FALSE
      }
    } else if (isTRUE(nrow(df) >= 3)) {
      # If there's 3 or more peps, keep those correlated
      # equal or better than median
      avgcorr <- rowMeans(df_cor)
      if (average.corr.cutoff == '>= median') {
        survivors <- names(avgcorr[avgcorr >= median(avgcorr)])
      }
      if (class(average.corr.cutoff) == 'numeric') {
        survivors <- names(avgcorr[avgcorr >= average.corr.cutoff])
      }
      indexpos <- data$peptide_id %in% survivors
      passed[indexpos] <- TRUE
    } else { # if only one peptide, do not keep
      passed[indexpos] <- FALSE
    }
  }
  result <- list(traces.wide=data[passed,], ids=traces$ids[passed,])
  class(result) <- 'traces'
  result #return the filtered traces data
}
>>>>>>> f46603dc81ae4a9ffbe532e5e5894bfa783ba38c
