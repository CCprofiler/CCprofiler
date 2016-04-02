#' Plot.
#' 
#' @export
plot.FragmentTraces <- function(trace.obj, plot=TRUE) {
    traces.long <- toLongFormat(trace.obj$traces)
    pl <- ggplot(traces.long) +
        ggtitle(paste(trace.obj$fragment_type, 'traces')) +
        xlab('fraction') + ylab('intensity')
    if (!hasAnalyteMap(trace.obj)) {
        pl <- pl + geom_line(aes(x=fraction, y=intensity, color=subunit.id.col))
    } else {
        pl <- pl + geom_line(aes(x=fraction, y=intensity, color=fragment_id,
                                 linetype=analyte_id))
    }
    if (plot) plot(pl)
    pl
}
