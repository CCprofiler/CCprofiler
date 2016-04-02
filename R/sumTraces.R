#' Produce analyte traces by summing traces of fragments. 
#' @param trace.obj An object of class \code{FragmentTraces}.
#' @return An object of type \code{FragmentTraces} belonging to a higher
#'         level of the compositional hierarchy.
#' @export
sumTraces <- function(trace.obj) {
    traces.long <- toLongFormat(trace.obj$traces)
    analyte.map <- trace.obj$analyte_map[, list(fragment_id, analyte_id)]
    traces.long.with.analyte <-
        merge(traces.long, analyte.map, by='fragment_id',
              allow.cartesian=TRUE)
    analyte.traces.long <-
        traces.long.with.analyte[, list(intensity=sum(intensity)),
                                 by=list(analyte_id, fraction)]
    setnames(analyte.traces.long, 'analyte_id', 'fragment_id')
    analyte.traces.wide <- toWideFormat(analyte.traces.long)
    trace.obj.copy <- trace.obj

    list(traces=analyte.traces.wide,
         fragment_type=trace.obj$analyte_type,
         analyte_type=NULL,
         analyte_map=NULL)
}
