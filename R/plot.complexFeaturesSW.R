#' Plot the result of the slidingWindow algorithm.
#' @export
plot.complexFeaturesSW <- function(sw.result, protein.traces.long, n.largest=NULL,
                          log=FALSE) {
    subgroups.dt <- sw.result$subgroups.dt
    complex.id <- unique(protein.traces.long$complex_id)
    trace.plot <- plotTraces(protein.traces.long,
                             'protein_id', 'complex_id',
                             paste('Protein traces of complex', complex.id),
                             plot=F,
                             log=log)
    if (!is.null(n.largest)) {
        n.subunits <- unique(subgroups.dt$n_subunits)
        n.subunits.ord <- n.subunits[order(n.subunits, decreasing=TRUE)]
        largest.n <- head(n.subunits.ord, n.largest)
        subgroups.dt <- subgroups.dt[n_subunits %in% largest.n]
    }

    subgroup.plot <- ggplot(subgroups.dt) +
        geom_point(aes(x=sec, y=protein_id, color=protein_id),
                  size=2, shape=15) +
        xlim(c(min(protein.traces.long$sec), max(protein.traces.long$sec))) +
        facet_wrap(~ subgroup, nrow=length(unique(subgroups.dt$subgroup))) +
        ggtitle(sprintf('Detected subgroups (window size = %d, correlation cutoff = %f)',
                         sw.result$window.size, sw.result$corr.cutoff)) +
        theme(axis.line.y=element_blank()
              # axis.text.y=element_blank(),
              # axis.ticks.y=element_blank(),
              # axis.title.y=element_blank()
              # panel.background=element_blank(),
              # panel.border=element_blank(),
              # panel.grid.major=element_blank(),
              # panel.grid.minor=element_blank(),
              # plot.background=element_blank()
              )
    p <- multiplot(trace.plot + ylim(c(0, 5e06)),
                   subgroup.plot)
    print(p)
    p
}

