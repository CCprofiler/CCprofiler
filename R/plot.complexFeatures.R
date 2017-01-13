#' Plot the result of the sliding window algorithm.
#' @import data.table
#' @import grid
#' @import gridExtra
#' @import gtable
#' @param sw.result An object of type 'complexFeaturesSW'.
#' @param trace.mat The same matrix that was passend to
#'        \code{findComplexFeatures}.
#' @param protein.names The names of proteins whose traces are contained in
#'        the matrix \code{trace.mat}.
#' @param n.largest Only plot the n largest subgroups. Optional.
#' @export
plot.complexFeatures <- function(res,
                                 traces.obj,
                                 complexID,
                                 calibration,
                                 pp=FALSE,
                                 sw=FALSE,
                                 estimate=FALSE,
                                 monomerMW=FALSE,
                                 log=FALSE) {


    features <- subset(res, complex_id == complexID)
    proteins <- unique(unlist(strsplit(features$subunits_annotated, split = ";")))
    traces <- subset(traces.obj,trace_ids = proteins)
    traces.long <- toLongFormat(traces$traces)

    detectedSubunits = strsplit(features$subunits_detected, ';')[[1]]
    sel_inComplex <- which(traces.long$id %in% detectedSubunits)
    traces.long$inComplex <- "FALSE"
    traces.long$inComplex[sel_inComplex] <- "TRUE"
    traces.long$inComplex = factor(traces.long$inComplex,levels=c("TRUE","FALSE"))
    # add monomer MW
    subunitMW = as.numeric(strsplit(features$monomer_sec, ';')[[1]])
    subunitMW.dt = data.table(id=detectedSubunits,mw=subunitMW)

    traces.long = merge(traces.long,subunitMW.dt,by="id",all.x = TRUE)

    up_map=uniprot_extended_table_human9606_20161130
    for (i in seq(1,nrow(traces.long),1)) {
      sel = which(up_map$Entry == traces.long$id[i])
      gene_names=unique(up_map$`Gene names  (primary )`[sel])
      if(length(gene_names) > 0) {
        traces.long$id[i] <- gene_names[1]
      }
    }

    setkey(traces.long, inComplex, id)

    complexName = unique(features$complex_name)[1]
    n_annotatedSubunits = features$n_subunits_annotated[1]
    n_signalSubunits = features$n_subunits_with_signal[1]
    n_detectedSubunits = features$n_subunits_detected[1]
    complexCompleteness = round(features$completeness[1],digits=2)

    grid.newpage()

    p <- ggplot(traces.long) +
      geom_line(aes_string(x='fraction', y='intensity', color='id',linetype='inComplex')) +
      ggtitle(paste0(complexName,
                     "\nAnnotated subunits: ",n_annotatedSubunits,
                     "   Subunits with signal: ",n_signalSubunits,
                     "\nCoeluting subunits: ",n_detectedSubunits,
                     "   Completeness: ",complexCompleteness)) +
      xlab('fraction') +
      ylab('intensity') +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()) +
      guides(fill=FALSE,linetype=FALSE) +
      theme(plot.margin = unit(c(1.5,.5,.5,.5),"cm")) +
      theme(plot.title = element_text(vjust=19,size=10))

    if (log) {
        p <- p + scale_y_log10('log(intensity)')
    }

#    p + geom_vline(aes(xintercept=left_sec, color=subgroup),
#                   data=found.features) +
#        geom_vline(aes(xintercept=right_sec, color=subgroup),
#                   data=found.features) +
#        theme(legend.position='none')

    if(pp==TRUE){
      p <- p + geom_rect(data=features,aes(xmin = left_pp, xmax = right_pp, ymin = -Inf, ymax = Inf, fill=subunits_detected),alpha = 0.25)
    }
    if(sw==TRUE){
      p <- p + geom_vline(data=features,aes(xintercept=left_sw), colour="grey",linetype="longdash")
      p <- p + geom_vline(data=features,aes(xintercept=right_sw), colour="grey",linetype="longdash")
    }
    if(estimate==TRUE){
      p <- p + geom_vline(data=features,aes(xintercept=complex_sec_estimated), colour="black",linetype="longdash")
    }
    if (monomerMW==TRUE){
      #p <- p + geom_vline(data=subset(traces.long,!is.na(mw)),aes(xintercept=mw,colour=id),linetype="solid")
      p <- p + geom_point(data = subset(traces.long,!is.na(mw)), mapping = aes(x = mw, y = Inf, colour=id),shape=18,size=5,alpha=.5)
      #p <- p + geom_segment(data=subset(traces.long,!is.na(mw)), aes(x=mw, y=0, xend=mw, yend=-Inf, colour=id))
    }
    p <- p + geom_vline(data=features,aes(xintercept=apex), colour="black",linetype="solid")


    p <- p + scale_x_continuous(name="SEC fraction",limits=c(0,80),breaks=seq(0,80,10),labels=seq(0,80,10))

    p2 <- p
    p2 <- p2 + scale_x_continuous(name="MW",
                   ,limits=c(0,80),breaks=seq(0,80,10),
                   labels=round(calibration$SECfractionToMW(seq(0,80,10)),digits=0)
                   )

    ## extract gtable
    g1 <- ggplot_gtable(ggplot_build(p))
    g2 <- ggplot_gtable(ggplot_build(p2))

    ## overlap the panel of the 2nd plot on that of the 1st plot
    pp <- c(subset(g1$layout, name=="panel", se=t:r))

    g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name=="panel")]], pp$t, pp$l, pp$b, pp$l)

    ## steal axis from second plot and modify
    ia <- which(g2$layout$name == "axis-b")
    ga <- g2$grobs[[ia]]
    ax <- ga$children[[2]]

    ## switch position of ticks and labels
    ax$heights <- rev(ax$heights)
    ax$grobs <- rev(ax$grobs)

    ## modify existing row to be tall enough for axis
    g$heights[[2]] <- g$heights[g2$layout[ia,]$t]

    ## add new axis
    g <- gtable_add_grob(g, ax, 2, 4, 2, 4)

    ## add new row for upper axis label
    g <- gtable_add_rows(g, g2$heights[1], 1)

    ## steal axis label from second plot
    ia2 <- which(g2$layout$name == "xlab-b")
    ga2 <- g2$grobs[[ia2]]
    g <- gtable_add_grob(g,ga2, 3, 4, 2, 4)

    # draw it
    #pdf("test.pdf",height=4,width=7)
    grid.draw(g)
    #dev.off()

}
