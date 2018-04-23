#' Add summary statistics to the fractionAnnotation object
#' @param traces Object of class traces.
#' @return Object of class traces with annotation of fractions.
#' @examples
#' inputTraces <- examplePeptideTraces
#' exampleTracesAnnotatedFractions <- annotateFractions(inputTraces)
#' exampleTracesAnnotatedFractions$fraction_annotation
#' @export

annotateFractions <- function(traces){
  UseMethod("annotateFractions", traces)
}

#' @describeIn annotateMolecularWeight Annotate single traces object
#' @export

annotateFractions.traces <- function(traces){
  .tracesTest(traces)
  traces_res <- copy(traces)
  traces_res <- annotateFractionMissingValues(traces)
  traces_res <- annotateFractionIdCounts(traces_res)
  traces_res <- annotateFractionIntensity(traces_res)
  traces_res <- annotateResidualIntensity(traces_res)
  traces_res <- annotateResidualIdCounts(traces_res)
  .tracesTest(traces_res)
  return(traces_res)
}

#' @describeIn annotateMolecularWeight Annotate multiple traces objects
#' @export

annotateFractions.tracesList <- function(traces){
  .tracesListTest(traces)
  res <- lapply(traces, annotateFractions.traces)
  class(res) <- "tracesList"
  .tracesListTest(res)
  return(res)
}

#' Calculates and annotates the number of missing values in every fraction
annotateFractionMissingValues <- function(traces, bound_left = 2, bound_right = 2, consider_borders = T){
  traces_imp <- findMissingValues(traces = traces, bound_left = bound_left, bound_right = bound_right, consider_borders = consider_borders)
  nas <- apply(getIntensityMatrix(traces_imp), 2, function(column) length(which(is.na(column))))
  traces$fraction_annotation$missingValues <- nas
  return(traces)
}

#' Calculates and annotates the number of ids in every fraction
annotateFractionIdCounts <- function(traces){
  .tracesTest(traces = traces)
  ids <- apply(getIntensityMatrix(traces), 2, function(column) length(which(column>0)))
  traces$fraction_annotation$nrIds <- ids
  return(traces)
}

#' Calculates and annotates the total intensity in every fraction
annotateFractionIntensity <- function(traces){
  .tracesTest(traces = traces)
  ints <- colSums(getIntensityMatrix(traces))
  traces$fraction_annotation$intensitySum <- ints
  return(traces)
}

#' Calculate loess residuals
annotateResidualIntensity <- function(traces){
  .tracesTest(traces = traces)
  ints <- colSums(getIntensityMatrix(traces))
  intSum <- data.frame(fraction = 1:length(ints), intensity = ints)
  l <- loess(intensity ~ fraction, data = intSum, span = 0.15)
  residuals <- residuals(l)
  traces$fraction_annotation$loessResiduals <- residuals
  return(traces)
}

#' Calculate loess residuals
annotateResidualIdCounts <- function(traces){
  .tracesTest(traces = traces)
  ints <- apply(getIntensityMatrix(traces), 2, function(column) length(which(column>0)))
  intSum <- data.frame(fraction = 1:length(ints), intensity = ints)
  l <- loess(intensity ~ fraction, data = intSum, span = 0.15)
  residuals <- residuals(l)
  traces$fraction_annotation$loessResidualsIdCounts <- residuals
  return(traces)
}
