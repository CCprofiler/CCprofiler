#' Get the shift between 2 traces objects via global corellation
#' @param traces1 First traces object to align.
#' @param traces2 Second traces object to align.
#' @param min_lag numeric, the lower bound of the range of lags to evaluate.
#' @param max_lag numeric, the upper bound of the range of lags to evaluate.
#' @param plot Logical, wether to plot a boxplot of corellation at different lags.
#' @return Numeric, the shift that produces the highest corellation
#' @export
getLag <- function(traces1, traces2, min_lag = -5, max_lag = 5, plot = T, ...){
  
  .tracesTest(traces1)
  .tracesTest(traces2)
  if(min_lag>max_lag) stop("the minimum lag must be smaller than the maximum lag!")
  if(min_lag>0 ) warning(paste0("Min lag is larger than zero, this will neccesarily result in a shift of at least ", min_lag))
  if( max_lag <0) warning(paste0("Max lag is smaller than zero, this will neccesarily result in a negative shift of at least ", min_lag))
  if(ncol(traces1$traces) != ncol(traces2$traces)){
    message("Traces objects have different number of fractions. Cropping the larger traces object for correlation calculation.")
    fracIntersect <- intersect(traces1$fraction_annotation$id, traces2$fraction_annotation$id)
  }else{
    fracIntersect <- NULL
  }
  
  ids1 <- traces1$traces[rowSums(traces1$traces[,!"id", with = F]) > 0, id]
  ids2 <- traces2$traces[rowSums(traces2$traces[,!"id", with = F]) > 0, id]
  common_ids <- intersect(ids1, ids2)
  
  traces1 <- subset(traces1, trace_subset_ids = common_ids, fraction_ids = fracIntersect)
  traces2 <- subset(traces2, trace_subset_ids = common_ids, fraction_ids = fracIntersect)
  setkey(traces1$traces, "id")
  setkey(traces2$traces, "id")
  # The order of the ids must be identical now -> numerical referencing is possible
  tr1 <- as.matrix(traces1$traces[,!"id", with = F])
  tr2 <- as.matrix(traces2$traces[,!"id", with = F])
  lag_range <- min_lag:max_lag    
  n <- ncol(tr1)
  cross_corr <- sapply(lag_range, function(l){
    t1 <- tr1[,max(1, l+1):min(n, n+l)]
    t2 <- tr2[,max(1, 1-l):min(n, n-l)]
    # sapply(1:nrow(t1), function(i) cor(t1[i,], t2[i,]))
    cA <- t1 - rowMeans(t1)
    cB <- t2 - rowMeans(t2)
    sA <- sqrt(rowMeans(cA^2))
    sB <- sqrt(rowMeans(cB^2))
    
    rowMeans(cA * cB) / (sA * sB)
  })
  # cross_corr <- cross_corr[-which(is.na(cross_corr), arr.ind = T)[,1],]
  cross_corr <- cross_corr[-which(is.nan(cross_corr), arr.ind = T)[,1],]
  if(plot){
    
    boxplot(cross_corr, names = lag_range, notch = T,
            xlab = "Fraction shift", ylab = "Correlation")
    lines(1:length(lag_range),colMeans(cross_corr))
    legend(x = length(lag_range) + 0.7, y = 1, legend = "mean", lty = 1, xjust = 1)
  }
  if(any(names(list(...)) == "returnPlot")){
    return(list(lag = lag_range[which.max(colSums(cross_corr))], 
                lag_range = lag_range,
                cross_corr = cross_corr))
  }else{
    return(lag_range[which.max(colSums(cross_corr))])
  }
  
  
}

#' Get the alignment shifts in tracesList object
#' @param tracesList Object of class tracesList to align.
#' @param alignment_order Integer vector, determines the sequence of alignment.
#' @param plot Logical, wether to plot boxplots of corellation at different lags.
#' @param PDF Logical, wether to produce a PDF in the working directory.
#' @param min_lag numeric, the lower bound of the range of lags to evaluate.
#' @param max_lag numeric, the upper bound of the range of lags to evaluate.
#' @param name Character string, the name of the PDF file.
#' @return Numeric vector, the shifts with the highest corellations for each alignment
#' @export

alignTraces <- function(tracesList,
                        alignment_order = seq_along(tracesList),
                        plot = TRUE,
                        PDF = FALSE,
                        min_lag = -5, max_lag = 5,
                        name= "AlignmentSummary.pdf", ...){
  .tracesListTest(tracesList)
  message("Sequentially aligning traces: ")
  if(PDF) pdf(gsub("$|\\.pdf$", ".pdf",name))
  lags <- sapply(alignment_order[-length(alignment_order)], function(i){
    message(paste(names(tracesList)[i], "to", names(tracesList)[i+1]))
    if(PDF|plot){
      l <- getLag(tracesList[[i]], tracesList[[i+1]], returnPlot = TRUE, plot = F,
                  min_lag = min_lag, max_lag = max_lag, ...)
      boxplot(l$cross_corr, names = l$lag_range, notch = T,
              xlab = "Fraction shift", ylab = "Correlation")
      title(paste(names(tracesList)[i], "to", names(tracesList)[i+1]))
      lines(1:length(l$lag_range),colMeans(l$cross_corr))
      legend(x = length(l$lag_range) + 0.7, y = 1, legend = "mean", lty = 1, xjust = 1)
      l$lag
    }else{
      getLag(tracesList[[i]], tracesList[[i+1]], plot = F,
             min_lag = min_lag, max_lag = max_lag, ...)
    }
  })
  if(PDF) dev.off()
  
  lags
}

#' Shift a traces object by n fractions
#' @param traces Object of class traces to shift
#' @param shift Numeric integer, how many fractions to shift. 
#' (positive = right shift, negative = left shift).
#' @param cropOverhang Logical, wether to remove the fraction(s) that is overhanging after shifting.
#' @param imputeMissing Logical, wether to impute the fraction(s) that is lagging after shifting.
#' @return traces object that is shifted by the specified fractions.
#' @export

shiftTraces <- function(traces, shift, cropOverhang = T, imputeMissing = T){
  .tracesTest(traces)
  tr <- copy(traces$traces)
  ids <- tr$id
  tr[, id := NULL]
  ncols <- ncol(tr)
  if(shift>0){
    for(i in 1:abs(shift)){
      trShift <- cbind(0,tr)
      if(imputeMissing){
        trShift[,1] <- pmax(0, 2 * trShift[[2]] - trShift[[3]])
      }
      trShift <- trShift[,1:ncols, with = F]
    }
    
  }else if(shift <0){
    for(i in 1:abs(shift)){
      trShift <- cbind(tr, 0)
      if(imputeMissing){
        trShift[,(ncols +1)] <- pmax(0, 2 * trShift[[ncols - 1]] - trShift[[ncols - 2]])
      }
      trShift <- trShift[,2:(ncols + 1), with = F]
    }
  }else if(shift == 0){
    return(traces)
  }else{
    stop("Shift must be numeric!")
  }
  
  names(trShift) <- as.character(1:ncols)
  trShift[,id := ids]
  traces$traces <- trShift
  .tracesTest(traces)
  return(traces)
}