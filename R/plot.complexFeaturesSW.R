#' Plot the result of the sliding window algorithm.
#' @param sw.result An object of type 'complexFeaturesSW'.
#' @param trace.mat The same matrix that was passend to
#'        \code{findComplexFeatures}.
#' @param protein.names The names of proteins whose traces are contained in
#'        the matrix \code{trace.mat}.
#' @param n.largest Only plot the n largest subgroups. Optional.
#' @export
plot.complexFeaturesSW <- function(sw.result,
                                   trace.mat,
                                   protein.names,
                                   n.largest=NULL,
                                   log=FALSE) {

    #@TODO change input structure!!!

    found.features <- sw.result$features

    # Produce a long list version of the trace matrix since its more convenient
    # for plotting.
    traces.dt <- data.table(trace.mat)
    traces.dt[, protein_id := protein.names]
    traces.long <- melt(traces.dt, id.vars='protein_id', value.name='intensity',
                        variable.name='fraction', variable.factor=FALSE)
    traces.long[, fraction := as.numeric(fraction)]
    setkey(traces.long, protein_id)


    p <- ggplot(traces.long) +
                geom_line(aes_string(x='fraction', y='intensity', color='protein_id')) +
                ggtitle('Protein traces') +
                xlab('fraction') +
                ylab('intensity')
    if (log) {
        p <- p + scale_y_log10('log(intensity)')
    }

    if (!is.null(n.largest)) {
        n.subunits <- unique(found.features$n_subunits)
        n.subunits.ord <- n.subunits[order(n.subunits, decreasing=TRUE)]
        largest.n <- head(n.subunits.ord, n.largest)
        found.features <- found.features[n_subunits %in% largest.n]
    }



#    p + geom_vline(aes(xintercept=left_sec, color=subgroup),
#                   data=found.features) +
#        geom_vline(aes(xintercept=right_sec, color=subgroup),
#                   data=found.features) +
#        theme(legend.position='none')

    p <- p + geom_rect(data=found.features,aes(xmin = left_sec, xmax = right_sec, ymin = -Inf, ymax = Inf,fill = subgroup),alpha = 0.25)
    p <- p + geom_vline(data=found.features,aes(xintercept=(log(mw_apparent) - 9.682387)/(-0.1043329), linetype = subgroup))
    #p <- p + geom_vline(data=found.features,aes(xintercept=(log(mw_estimated) - 9.682387)/(-0.1043329), linetype = subgroup))

}
