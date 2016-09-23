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
