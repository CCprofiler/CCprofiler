#' annotateMassDistribution
#' @description Estimates fraction of assembled vs. monomeric mass per id.
#' @param traces An object of type traces. These should be the raw, unprocessed traces.
#' @export
annotateMassDistribution <- function(traces){
  UseMethod("annotateMassDistribution", traces)
}

#' @describeIn annotateMassDistribution Estimates fraction of assembled vs. monomeric mass per id.
#' @export
annotateMassDistribution.traces <- function(traces){
  .tracesTest(traces)

  if(! ("protein_mw" %in% names(traces$trace_annotation))) {
    stop("annotateMassDistribution requires information on the monomer molecular weight. Please see annotateTraces.")
  }

  if(! ("molecular_weight" %in% names(traces$fraction_annotation))) {
    stop("annotateMassDistribution requires information on the molecular weight assigned to each fraction. Please see annotateMolecularWeight.")
  }

  intMat <- getIntensityMatrix(traces)
  mw_info <- traces$trace_annotation
  max_trace <- max(traces$fraction_annotation$id)
  mw_info[,assembly_boundary := 2*protein_mw]
  mw_info[is.na(assembly_boundary), assembly_boundary := 0]
  mw_info[,assembly_boundary_fraction := unlist(lapply(assembly_boundary,function(x){max(traces$fraction_annotation[molecular_weight>=x]$id)}))]

  # Get intensity sum
  mw_info[,sum_assembled := sum(intMat[id,1:assembly_boundary_fraction-1]), id]
  mw_info[,sum_monomeric := sum(intMat[id,assembly_boundary_fraction:max_trace]), id]
  mw_info[,sum := sum(intMat[id,]), id]
  mw_info[is.na(assembly_boundary), sum_assembled := 0]
  mw_info[is.na(assembly_boundary), sum_monomeric := 0]
  mw_info[, sum_assembled_norm := 1/sum*sum_assembled]
  mw_info[, sum_monomeric_norm := 1/sum*sum_monomeric]
  #mw_info[, test := sum_monomeric_norm+sum_assembled_norm]

  traces$trace_annotation <- mw_info
  .tracesTest(traces)
  return(traces)
}

#' @describeIn annotateMassDistribution Estimates fraction of assembled vs. monomeric mass per id.
#' @export
annotateMassDistribution.tracesList <- function(tracesList){
  .tracesListTest(tracesList)
  tracesListRes <- lapply(tracesList, annotateMassDistribution.traces)
  class(tracesListRes) <- "tracesList"
  .tracesListTest(tracesListRes)
  return(tracesListRes)
}

getMassAssemblyChange <- function(tracesList){
  .tracesListTest(tracesList)

}
