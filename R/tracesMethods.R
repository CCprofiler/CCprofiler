# Due to: http://stackoverflow.com/questions/24501245/data-table-throws-object-not-found-error
.datatable.aware=TRUE

#' Subset traces by trace_ids or fraction_ids
#' @description  Subset a taces object by trace ids and/or fraction ids.
#' @param traces Object of class traces.
#' @param trace_ids Character vector specifying the trace identifiers
#'        for subsetting traces, e.g. peptide or protein ids.
#' @param fraction_ids Numeric vector specifying the fraction identifiers
#'        for subsetting traces. The resulting traces object will not have
#'        fraction ids starting from 1, co-elution feature finding might thus
#'        be compromised. Please don't use this option for general
#'        pre-processing purposes.
#' @return traces Object of class traces.
#' @examples
#' # Load some example data:
#' inputTraces <- examplePeptideTraces
#' subsetPeptides <- c("AIIDEFEQK","AIQLSGAEQLEALK","AKEALIAASETLK")
#' subsetFractions <- c(5:20)
#' # Run subsetting:
#' subsettedTraces <- subset(inputTraces,trace_ids=subsetPeptides,fraction_ids=subsetFractions)
#' # Inspect subsetting result:
#' subsettedTraces
#' summary(subsettedTraces)
#' @export
subset.traces <- function(traces,trace_ids=NULL,fraction_ids=NULL){
    if (!is.null(trace_ids)) {
      traces$traces <- subset(traces$traces,id %in% trace_ids)
      traces$trace_annotation <- subset(traces$trace_annotation, id %in% trace_ids)
    }
    if (!is.null(fraction_ids)){
      traces$traces <- subset(traces$traces,select=c(names(traces$traces)[fraction_ids],"id"))
      traces$fraction_annotation <- subset(traces$fraction_annotation ,id %in% fraction_ids)
    }
    traces
}

#' Get intensity matrix from traces object
#' @description  Get a matrix of intensity values from a traces object.
#' @param traces Object of class traces.
#' @return Numeric matrix with intensity values.
#' @examples
#' intensityMatrix <- getIntensityMatrix(examplePeptideTraces)
#' head(intensityMatrix)
#' @export
getIntensityMatrix <- function(traces) {
    ids <- traces$traces$id
    intensity.mat <- as.matrix(sapply(subset(traces$traces,
                                      select=-id),as.numeric))
    rownames(intensity.mat) <- ids
    intensity.mat
}


#' Convert traces from wide format to long format
#' @description Convert the data.table traces from a trace object from wide format to long format.
#' @param traces.dt A data.table with an id column \code{id} and
#'        columns of continuously numbered fractions.
#' @return A data.table with columns
#'          \itemize{
#'           \item \code{id}
#'           \item \code{fraction}
#'           \item \code{intensity}
#'          }
#' @examples
#' tracesWide <- examplePeptideTraces$traces
#' tracesLong <- toLongFormat(tracesWide)
#' tracesLong
#' @export
toLongFormat <- function(traces.dt) {
  traces.dt.long <-
    melt(traces.dt, id.var='id', variable.name='fraction',
         value.name='intensity', variable.factor=FALSE)
  traces.dt.long[, fraction := as.numeric(fraction)]
  setkey(traces.dt.long,id)
  data.table(traces.dt.long)
  traces.dt.long
}

#' Plot traces
#' @description Plot a traces object as chromatograms.
#' @param traces Object of class traces.
#' @param legend Logical, whether a legend of the traces should be plotted. Default is \code{TRUE}.
#' @param PDF Logical, whether to plot to PDF. Default is \code{FALSE}.
#' @param name Character string with name of the plot, only used if \code{PDF=TRUE}. Default is "Traces".
#' @examples
#' plot(exampleProteinTraces)
#' @export
plot.traces <- function(traces, ledgend = TRUE, PDF=FALSE, name="Traces") {
  traces.long <- toLongFormat(traces$traces)
  pl <- ggplot(traces.long)
  pl <- pl + ggtitle(paste(traces$trace_type, 'traces'))
  pl <- pl + xlab('fraction') + ylab('intensity')
  pl <- pl + geom_line(aes(x=fraction, y=intensity, color=id))
  if (!ledgend) {
    pl <- pl + theme(legend.position="none")
  }
  if(PDF){
    pdf(paste0(name,".pdf"))
  }
  plot(pl)
  if(PDF){
    dev.off()
  }
}


#' Summarize a traces object
#' @description Summarize a traces object to get an overview.
#' @param traces Object of class traces.
#' @return A summary list with following entries:
#'          \itemize{
#'           \item \code{metrics} Summary of traces counts including decoy statistics.
#'           \item \code{type} Type of traces: peptides or proteins.
#'           \item \code{annotations} Names of trace annotations.
#'           \item \code{fraction_count} Number of fractions in traces object.
#'           \item \code{SibPepCorr} Summary of sibling peptide correlation.
#'           Only reported if sibling peptide correlations were previously calculated.
#'          }
#' @examples
#' summary(examplePeptideTracesFiltered)
#' @export
summary.traces <- function(traces) {
  no_traces <- nrow(traces$trace_annotation)
  no_decoys <- length(grep("DECOY", traces$trace_annotation$protein_id))
  no_targets <- no_traces - no_decoys
  pct_decoys <- signif(no_decoys/no_traces * 100, 2)
  res <- c(no_traces, no_targets, no_decoys, pct_decoys)
  names(res) <- c("No. of Traces", "No. of Targets", "No. of Decoys", "% Decoys")
  annotation_info <- names(traces$trace_annotation)
  fraction_info = length(traces$fraction_annotation$id)
  type=traces$trace_type
  if("SibPepCorr" %in% names(traces$trace_annotation)) {
    SibPepCorr_summary <- summary(traces$trace_annotation$SibPepCorr)
    summary=list(metrics=res,type=type,annotations=annotation_info,fraction_count=fraction_info,SibPepCorr_summary=SibPepCorr_summary)
  } else {
    summary=list(metrics=res,type=type,annotations=annotation_info,fraction_count=fraction_info)
  }
  summary
}
