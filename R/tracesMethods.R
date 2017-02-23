#' Subset a traces.obj by trace_ids or fraction_ids.
#' @param traces.obj An object of type \code{traces.obj}.
#' @param trace_ids A character vector specifying the trace identifiers
#'        for subsetting \code{traces.obj}.
#' @param fraction_ids A numeric vector specifying the fraction identifiers
#'        for subsetting \code{traces.obj}.
#' @return traces.obj An object of type \code{traces.obj}.
#' @export
subset.traces <- function(traces.obj,trace_ids=NULL,fraction_ids=NULL){
    if (!is.null(trace_ids)) {
      traces.obj$traces <- subset(traces.obj$traces,id %in% trace_ids)
      traces.obj$trace_annotation <- subset(traces.obj$trace_annotation ,id %in% trace_ids)
    }
    if (!is.null(fraction_ids)){
      traces.obj$traces <- traces.obj$traces[,c(fraction_ids,"id"),with=FALSE]
      traces.obj$fraction_annotation <- subset(traces.obj$fraction_annotation ,id %in% fraction_ids)
    }
    traces.obj
}

#' Get a matrix of intensity values from a traces object.
#' @param traces.obj An object of type \code{traces.obj}.
#' @return A matrix with intensity values.
#' @export
getIntensityMatrix <- function(traces.obj) {
    ids <- traces.obj$traces$id
    intensity.mat <- as.matrix(sapply(subset(traces.obj$traces,
                                      select=-id),as.numeric))
    rownames(intensity.mat) <- ids
    intensity.mat
}


#' Convert a data.table containing traces from wide format to long format.
#' @param traces.dt A data.table with an id column \code{id} and
#'        columns of continuously numbered fractions.
#' @return A data.table with columns
#'          \itemize{
#'           \item \code{id}
#'           \item \code{fraction}
#'           \item \code{intensity}
#'          }
toLongFormat <- function(traces.dt) {
  traces.dt.long <-
    melt(traces.dt, id.var='id', variable.name='fraction',
         value.name='intensity', variable.factor=FALSE)
  traces.dt.long[, fraction := as.numeric(fraction)]
  setkey(traces.dt.long,id)
  data.table(traces.dt.long)
  traces.dt.long
}

#' Plot a traces.obj.
#' @param traces.obj An object of type \code{traces.obj}.
#' @export
plot.traces <- function(traces.obj, plot=TRUE, ledgend = TRUE, title=NULL) {
  traces.long <- toLongFormat(traces.obj$traces)
  pl <- ggplot(traces.long)
  if (is.null(title)) {
    pl <- pl + ggtitle(paste(traces.obj$trace_type, 'traces'))
  } else {
    pl <- pl + ggtitle(title)
  }
  pl <- pl + xlab('fraction') + ylab('intensity')
  pl <- pl + geom_line(aes(x=fraction, y=intensity, color=id))
  if (!ledgend) {
    pl <- pl + theme(legend.position="none")
  }
  if (plot) plot(pl)
  pl
}


#' Summarize a \code{traces.obj}
#' @param traces.obj An object of type \code{traces.obj}.
#' @export
summary.traces <- function(traces.obj) {
  no_traces <- nrow(traces.obj$trace_annotation)
  no_decoys <- length(grep("DECOY", traces.obj$trace_annotation$protein_id))
  no_targets <- no_traces - no_decoys
  pct_decoys <- signif(no_decoys/no_traces * 100, 2)
  res <- c(no_traces, no_targets, no_decoys, pct_decoys)
  names(res) <- c("No. of Traces", "No. of Targets", "No. of Decoys", "% Decoys")
  res
}

