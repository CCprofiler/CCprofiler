#' Plot protein features.
#'
#' @param res proteinFeatures
#' @param peptideTraces An object of type \code{traces}.
#' @param proteinID protein ID
#' @param calibration list of calibraion functions form calibrateSECMW
#' @param window.size Size of the window. Numeric.
#' @param plot_peak logical
#' @param plot_monomer logical
#' @export

plotProteinFeatures <- function(res,
                                 peptideTraces,
                                 proteinID,
                                 calibration,
                                 plot_peak=TRUE,
                                 plot_monomer=TRUE,
                                 log=FALSE,
                                 plot_apex=TRUE,
                                 plot_in_complex_estimate=TRUE) {


  features <- subset(res, protein_id == proteinID)
  #peptides <- unique(unlist(strsplit(features$subunits_annotated, split = ";")))
  peptides <- peptideTraces$trace_annotation[protein_id == proteinID]$id
  traces <- subset(peptideTraces,trace_ids = peptides)
  traces.long <- toLongFormat(traces$traces)

  # add monomer MW
  monomer = as.numeric(unlist(features$monomer_sec))
  monomer_max = as.numeric(unlist(features$monomer_sec_max))

  setkey(traces.long, id)

  proteinName = unique(features$protein_id)[1]

  group.colors <- c("TRUE" = "green3", "FALSE" = "firebrick1")

  grid.newpage()

  p <- ggplot(traces.long) +
    geom_line(aes_string(x='fraction', y='intensity', color='id')) +
    ggtitle(proteinName) +
    xlab('fraction') +
    ylab('intensity') +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
    theme(plot.margin = unit(c(1,.5,.5,.5),"cm")) +
    theme(plot.title = element_text(vjust=19,size=10)) +
    guides(color=FALSE)

    if (log) {
      p <- p + scale_y_log10('log(intensity)')
    }

  if(plot_peak==TRUE){
    if (plot_in_complex_estimate == TRUE) {
      p <- p + geom_rect(data=features,aes(xmin = left_pp, xmax = right_pp, ymin = 0, ymax = Inf, fill=in_complex),alpha = 0.25) +
      scale_fill_manual(values=group.colors)
    } else {
      p <- p + geom_rect(data=features,aes(xmin = left_pp, xmax = right_pp, ymin = 0, ymax = Inf),alpha = 0.25)
    }
    p <- p + geom_vline(data=features,aes(xintercept = left_pp), colour="darkgrey", linetype="dashed")
    p <- p + geom_vline(data=features,aes(xintercept = right_pp), colour="darkgrey",linetype="dashed")
  }

  if (plot_monomer==TRUE){
    p <- p + geom_vline(xintercept=monomer,linetype="solid")
    p <- p + geom_vline(xintercept=monomer_max,linetype="dashed")
  }
  if(plot_apex==TRUE){
    p <- p + geom_vline(data=features,aes(xintercept=apex), colour="darkgrey", linetype="solid",size=1.5)
  }
  #p <- p + theme(legend.position="none")
  #print(p)
  p2 <- p
  p <- p + scale_x_continuous(name="SEC fraction",limits=c(0,80),breaks=seq(0,80,10),labels=seq(0,80,10))
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
}
