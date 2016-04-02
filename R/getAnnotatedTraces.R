#' Get a data.table of intensities with two columns \code{analyte_id} and
#' \code{fragment_id}.
#' @param trace.obj An object of type \code{FragmentTraces}.
#' @return A data.table of annotated traces.
#' @export
getAnnotatedFragmentTraces <- function(trace.obj) {
    analyte.map <- trace.obj$analyte_map[, list(fragment_id, analyte_id)]
    merged.traces <- merge(trace.obj$traces, analyte.map, 
                           by='fragment_id', allow.cartesian=TRUE)
    column.names <- colnames(merged.traces)
    fraction.columns <-
        setdiff(column.names, c('analyte_id', 'fragment_id'))
    setcolorder(merged.traces, c(c('analyte_id', 'fragment_id'),
                                 fraction.columns))
    merged.traces
}
