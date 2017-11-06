#' Normalize trace intensities to reference traces
#' @param traces Object of class traces or tracesList.
#' @param standard_ids Character vector, peptide or protein ids of reference peptides
#' @param method Character string, whether to normalize to the median or the mean intensity 
#' of the standard.
#' @return Object of class traces with normalized intensities.
#' @export

normalizeToStandard <- function(traces, standard_ids, method = "median",
                                filter_complete = TRUE, top_n = NULL){
  UseMethod("normalizeToStandard", traces)
}

#' @describeIn normalizeToStandard Normalize to reference traces in multiple traces objects
#' @export

normalizeToStandard.traces <- function(traces, standard_ids, method = "median",
                                       filter_complete = TRUE, top_n = NULL, ...){
  # .tracesTest(traces)
  if(any(standard_ids %in% traces$trace_annotation$protein_id)){
    standards <- getIntensityMatrix(subset(traces,
                                           trace_subset_ids = standard_ids,
                                           trace_subset_type = "protein_id"))
  }else{
    standards <- getIntensityMatrix(subset(traces, trace_subset_ids = standard_ids))
  }
  if(filter_complete){
    completeRows <- !apply(standards, 1, function(x) any(x == 0))
    if(any(completeRows)){
      standards <- standards[completeRows, ]
    }else{
      stop("No complete observations of the standard found")
    }
    
  }
  if(!is.null(top_n)){
    standards <- standards[order(rowSums(standards), decreasing = T),]
    n <- nrow(standards)
    if(n < top_n) message(paste("Warning: Only", n, "Standards found. Using those..."))
    standards <- standards[1:min(top_n, n),]
  }
  if(method == "median"){
    fractionMedians <- apply(standards, 2, function(x) median(x[x>0]))
    medianOfMedians <- median(fractionMedians)
    normVal <- fractionMedians / medianOfMedians
    normTraces <- normalizeTraces(traces, normVal, method = "/")
  }else if(method == "mean"){
    fractionMeans <- apply(standards, 2, function(x) mean(x[x>0]))
    meanOfMeans <- mean(fractionMeans)
    normVal <- fractionMeans / meanOfMeans
    normTraces <- normalizeTraces(traces, normVal, method = "/")
  }else{
    stop(paste(method, "is not a valid method. Must be either 'median' or 'mean'."))
  }
  return(normTraces)
}

#' @describeIn normalizeToStandard Normalize to reference traces in multiple traces objects
#' @export

normalizeToStandard.tracesList <- function(tracesList, standard_ids, method = "median",
                                           filter_complete = TRUE, top_n = NULL, ...){
  .tracesListTest(tracesList)
  res <- lapply(tracesList, normalizeToStandard.traces,
                standard_ids = standard_ids, method = method,
                filter_complete = filter_complete, top_n = top_n, ...)
  class(res) <- "tracesList"
  .tracesListTest(res)
  return(res)
}


#' Normalize trace intensities by spline smoothing over the fractions
#' @description Fits a spline to a summary statistic of every elution fraction
#' and normalizes each fraction to fit the spline
#' @param method Character string, whether to normalize to the median or the mean intensity 
#' @param scale Character string wether to transform values before smoothing
#' @param span Numeric, controls the degree of smoothing (see also \code[loess] function)
#' @param plot Logical, whether to plot a Visualization of the smoothing
#' @param smoothe_on_targets Logical, whether to exclude decoy traces from smoothing
#' @param PDF Logical, whether to produce a PDF of the summary plot
#' @param name Character string, the name of the PDF produced. 
#' @return Object of class traces with normalized intensities.
#' @export


smootheTraces <- function(traces,
                          method = c("mean", "median"),
                          scale = c("none", "log"),
                          span = 0.15,
                          plot = FALSE,
                          PDF = F,
                          name = "Smoothing_summary.pdf",
                          smoothe_on_targets = TRUE, ...){
  UseMethod("smootheTraces", traces)
}

#' @describeIn smootheTraces Smoothe traces of a single traces object
#' @export
smootheTraces.traces <- function(traces,
                          method = c("mean", "median"),
                          scale = c("none", "log"),
                          span = 0.15,
                          plot = FALSE,
                          PDF = F,
                          name = "Smoothing_summary.pdf",
                          smoothe_on_targets = TRUE, ...){
  # .tracesTest(traces)
  method <- match.arg(method)
  scale <- match.arg(scale)
  
  ## Get intensity values
  targetIds <- traces$trace_annotation[!grepl("DECOY_", protein_id), id]
  if(smoothe_on_targets == TRUE){
    targetTraces <- subset(traces, trace_subset_ids = targetIds)
    targetInt <- getIntensityMatrix(targetTraces)
  }else{
    targetInt <- getIntensityMatrix(traces)
  }
  
  
  ## Scale
  smallVal <- targetInt != 0 & targetInt < 1
  if(any(smallVal)){
    targetInt[smallVal] <- 1
    message("There are intensities which are less than 1. These intensities are replaced with 1.")
  }
  
  if(scale == "log"){
    targetInt <- log2(targetInt)
  }else if (scale != "none"){
    stop("Invalid parameter 'scale' must be either 'log' or 'none")
  }
  
  ## Calculate fraction summaries
  if(method == "median"){
    ints <- apply(targetInt, 2, function(x) median(x[x > 0]))
    # }else if(method == "sum"){
    #   ints <- colSums(targetInt)
  }else if(method == "mean"){
    ints <- apply(targetInt, 2, function(x) mean(x[x > 0]))
  }
  
  intSum <- data.frame(fraction = 1:ncol(targetInt), intensity = ints)
  
  
  ## Fit a local regression to smooth the fractions 
  l <- loess(intensity ~ fraction, data = intSum, span = span)
  intSum$model <- predict(l)
  intSum$residuals <- residuals(l)
  
  ## Plot smoothing result
  p <- ggplot(intSum, aes(x = fraction, y = intensity)) +
    geom_point() +
    geom_line(aes(y = model)) +
    theme_bw()
  
  if(PDF) pdf(gsub("$|\\.pdf$", ".pdf",name))
  if(plot){
    plot(p)
  }
  if(PDF) dev.off()
  if(scale == "log"){
    normTraces <- normalizeTraces(traces, 2^intSum$residuals,
                                  method = "/")
    
  }else if(scale == "none"){
    normTraces <- normalizeTraces(traces, intSum$intensity/intSum$model, method = "/")
  }
  if(any(names(list(...)) == "returnPlot")){
    return(list(traces = normTraces, summaryPlot = p))
  }else{
    return(normTraces)
  }
  
}
#' @describeIn smootheTraces Smoothe traces of a tracesList object
#' @export

smootheTraces.tracesList <- function(traces,
                                 method = c("mean", "median"),
                                 scale = c("none", "log"),
                                 span = 0.15,
                                 plot = FALSE,
                                 PDF = F,
                                 name = "Smoothing_summary.pdf",
                                 smoothe_on_targets = TRUE, ...){
  .tracesListTest(traces)
  
  if(PDF) pdf(gsub("$|\\.pdf$", ".pdf",name))
  res <- lapply(names(traces), function(tr){
    r <- smootheTraces.traces(traces = traces[[tr]],                                
                              method = method,
                              scale = scale,
                              span = span,
                              plot =F,
                              PDF = F,
                              name = name,
                              smoothe_on_targets = smoothe_on_targets,
                              returnPlot = TRUE)
    p <- r[["summaryPlot"]] + ggtitle(paste0("Smoothing summary - ",tr))
    plot(p)
    return(r[["traces"]])
  })
  if(PDF) dev.off()
  names(res) <- names(traces)
  class(res) <- "tracesList"
  .tracesListTest(res)
  return(res)
}

  
normalizeTraces <- function(traces, normVals, method = "-", transform = c("none", "log", "2^")){
  transform <- match.arg(transform)
  normTraces <- copy(traces)
  normTraces$traces[,id := NULL]
  normTraces$traces <- as.data.table(sweep(normTraces$traces, 2, normVals, method))
  if(transform == "2^"){
    normTraces$traces <- as.data.table(2^normTraces$traces)
  }else if(transform == "log"){
    normTraces$traces <- as.data.table(log2(normTraces$traces))
  }
  # normTraces$traces[normTraces$traces < 0] <- 0
  normTraces$traces[,id := traces$traces$id]
  return(normTraces)
}



