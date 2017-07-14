#' Plot the result of the sliding window algorithm.
#' @description Chomatographic visualization
#' @import data.table
#' @import grid
#' @import gridExtra
#' @import gtable
#' @param res data.table
#' @param proteinTraces traces object of type protein
#' @param complexID character string specifying complex_id
#' @param annotationID character string specifying column name in trace annotation to use for trace labeling, default is protein_id (uniprot id)
#' @param peak_area logical if selected peak area should be highlighted
#' @param sliding_window_are logical if original sliding_window area should be highlighted
#' @param estimated_complex_MW logical if estimated complex MW should be indicated
#' @param monomer_MW logical if monomer MWs should be indicated
#' @param log logical if intensities should be log transformed
#' @export
plotComplexFeatures <- function(res,
                               traces,
                               complexID,
                               annotationID="protein_id",
                               apex=TRUE,
                               peak_area=FALSE,
                               sliding_window_area=FALSE,
                               estimated_complex_MW=FALSE,
                               monomer_MW=FALSE,
                               log=FALSE,
                               legend = TRUE,
                               PDF=FALSE) {
  .tracesTest(traces)
  features <- subset(res, complex_id == complexID)
  proteins <- unique(unlist(strsplit(features$subunits_annotated, split = ";")))
  traces <- subset(traces, trace_subset_ids = proteins)
  traces.long <- toLongFormat(traces$traces)
  traces.long <- merge(traces.long,traces$fraction_annotation,by.x="fraction",by.y="id")
  detectedSubunits = strsplit(features$subunits_detected, ';')[[1]]
  sel_inComplex <- which(traces.long$id %in% detectedSubunits)
  traces.long$inComplex <- "FALSE"
  traces.long$inComplex[sel_inComplex] <- "TRUE"
  traces.long$inComplex = factor(traces.long$inComplex,levels=c("TRUE","FALSE"))
  # add monomer MW
  subunitMW = as.numeric(strsplit(features$monomer_sec, ';')[[1]])
  subunitMW.dt = data.table(id=detectedSubunits,mw=subunitMW)
  
  complexName = unique(features$complex_name, by="inComplex")[1]
  n_annotatedSubunits = features$n_subunits_annotated[1]
  n_signalSubunits = features$n_subunits_with_signal[1]
  n_detectedSubunits = features$n_subunits_detected[1]
  complexCompleteness = round(features$completeness[1],digits=2)
  p <- ggplot(traces.long) +
    geom_line(aes_string(x='fraction', y='intensity', color='id')) +
    xlab('fraction') +
    ylab('intensity') +
    theme_bw() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5),"cm")) +
    ggtitle(paste0(complexName,
                   "\nAnnotated subunits: ",n_annotatedSubunits,
                   "   Subunits with signal: ",n_signalSubunits,
                   "\nCoeluting subunits: ",n_detectedSubunits,
                   "   Completeness: ",complexCompleteness)) +
    theme(plot.title = element_text(vjust=19,size=10, face="bold")) +
    guides(fill=FALSE)
  if (log) {
    p <- p + scale_y_log10('log(intensity)')
  }
  if(peak_area==TRUE){
    p <- p + geom_rect(data=features,aes(xmin = left_pp, xmax = right_pp, ymin = -Inf, ymax = Inf, fill=subunits_detected),alpha = 0.25)
  }
  if(sliding_window_area==TRUE){
    p <- p + geom_vline(data=features,aes(xintercept=left_sw), colour="grey",linetype="longdash")
    p <- p + geom_vline(data=features,aes(xintercept=right_sw), colour="grey",linetype="longdash")
  }
  if(estimated_complex_MW==TRUE){
    p <- p + geom_vline(data=features,aes(xintercept=complex_sec_estimated), colour="black",linetype="longdash")
  }
  if (monomer_MW==TRUE){
    p <- p + geom_point(data = subunitMW.dt, mapping = aes(x = subunitMW, y = Inf, colour=id),shape=18,size=5,alpha=.5)
  }
  if (apex==TRUE){
    p <- p + geom_vline(data=features,aes(xintercept=apex), colour="black",linetype="solid")
  }
  if (!legend) {
    p <- p + theme(legend.position="none")
  }
  
  if ("molecular_weight" %in% names(traces$fraction_annotation)) {
    p2 <- p
    p <- p + scale_x_continuous(name="fraction",
                                breaks=seq(min(traces$fraction_annotation$id),
                                           max(traces$fraction_annotation$id),10),
                                labels=seq(min(traces$fraction_annotation$id),
                                           max(traces$fraction_annotation$id),10))
    p2 <- p2 + scale_x_continuous(name="molecular weight (kDa)",
                                  breaks=seq(min(traces$fraction_annotation$id),max(traces$fraction_annotation$id),10),
                                  labels=round(traces$fraction_annotation$molecular_weight,digits=0)[seq(1,length(traces$fraction_annotation$id),10)]
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
    
    if(PDF){
      pdf(paste0(name,".pdf"))
    }
    grid.draw(g)
    if(PDF){
      dev.off()
    }
  }else{
    if(PDF){
      pdf(paste0(name,".pdf"))
    }
    plot(p)
    if(PDF){
      dev.off()
    }
  }
}
