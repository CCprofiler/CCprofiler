#' Get a data.table of annotated Traces
#' @param trace.obj An object of type \code{Traces}.
#' @return A data.table of annotated traces.
#' @export
getAnnotatedTraces <- function(traces.obj) {
    merged.traces <- merge(traces.obj$traces, traces.obj$trace_annotation, by='id', allow.cartesian=TRUE)
    setorder(merged.traces, "id")
    merged.traces
}
