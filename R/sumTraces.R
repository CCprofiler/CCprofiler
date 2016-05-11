#' Produce analyte traces by summing traces of fragments. 
#' @param trace.obj An object of class \code{FragmentTraces}.
#' @param by The annotation column that should be used when grouping fragments
#'        for summation.
#' @return An object of type \code{FragmentTraces} with summed traces.
#' @export
sumTraces <- function(trace.obj, by='protein_id') {
    # traces.long <- toLongFormat(trace.obj$traces)
    # analyte.map <- trace.obj$analyte_map[, list(fragment_id, analyte_id)]
    # traces.long.with.analyte <-
    #     merge(traces.long, analyte.map, by='fragment_id',
    #           allow.cartesian=TRUE)
    # analyte.traces.long <-
    #     traces.long.with.analyte[, list(intensity=sum(intensity)),
    #                              by=list(analyte_id, fraction)]
    # setnames(analyte.traces.long, 'analyte_id', 'fragment_id')
    # analyte.traces.wide <- toWideFormat(analyte.traces.long)
    # trace.obj.copy <- trace.obj

    # list(traces=analyte.traces.wide,
    #      fragment_type=trace.obj$analyte_type,
    #      analyte_type=NULL,
    #      analyte_map=NULL)
}
