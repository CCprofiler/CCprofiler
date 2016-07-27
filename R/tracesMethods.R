
#' Subset a traces.obj by trace_ids or fraction_ids.
#' @param traces.obj An object of type \code{traces.obj}.
#' @param trace_ids A character vector specifying the trace identifiers
#'        for subsetting \code{trace_ids}.
#' @param fraction_ids A numeric vector specifying the fraction identifiers
#'        for subsetting \code{fragment_ids}.
#' @return traces.obj An object of type \code{FragmentTraces}.
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

#' Get a matrix of fragment intensity values.
#' @param traces.obj An object of type \code{traces.obj}.
#' @return A matrix of intensity values.
#' @export
getIntensityMatrix.traces <- function(traces.obj) {
    ids <- traces.obj$traces$id
    intensity.mat <- as.matrix(subset(traces.obj$traces,
                                      select=-id))
    rownames(intensity.mat) <- ids
    intensity.mat
}


#' Convert a data.table containing traces from wide format to long format.
#' @param traces.dt A data.table with an id column \code{fragment_id} and
#'        columns of continuously numbered fractions.
#' @return A data.table with columns
#'          \itemize{
#'           \item \code{fragment_id}
#'           \item \code{fraction}
#'           \item \code{intensity}
#'          }
toLongFormat <- function(traces.dt) {
    traces.dt.long <-
        melt(traces.dt, id.var='id', variable.name='fraction',
             value.name='intensity', variable.factor=FALSE)
    traces.dt.long[, fraction := as.numeric(fraction)]
    print(traces.dt.long) #@TODO fix bug
    traces.dt.long
}

#' Plot.
#' @param traces.obj An object of type \code{traces.obj}.
#' @export
plot.traces <- function(traces.obj, plot=TRUE) {
    traces.long <- toLongFormat(traces.obj$traces)
    pl <- ggplot(traces.long) +
        ggtitle(paste(traces.obj$trace_type, 'traces')) +
        xlab('fraction') + ylab('intensity')
        pl <- pl + geom_line(aes(x=fraction, y=intensity, color=id))
    if (plot) plot(pl)
    pl
}
