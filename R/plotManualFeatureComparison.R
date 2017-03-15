#' plotManualFeatureComparison
#' @description plotManualFeatureComparison
#' @import data.table
#' @import grid
#' @import gridExtra
#' @import gtable
#' @param mapped_features proteinFeatures
#' @param traces An object of type \code{traces}.
#' @param manual_stats Stats of manual annotation.
#' @param calibration list of calibraion functions form calibrateSECMW
#' @param proteinID protein ID
#' @param plot_peak TRUE,
#' @param plot_monomer TRUE,
#' @param log FALSE,
#' @param plot_apex TRUE,
#' @param plot_in_complex_estimate TRUE,
#' @param TP FALSE
#' @export
plotManualFeatureComparison <- function(mapped_features,
                                 traces,
                                 manual_stats,
                                 calibration,
                                 proteinID,
                                 plot_peak=TRUE,
                                 plot_monomer=TRUE,
                                 log=FALSE,
                                 plot_apex=TRUE,
                                 plot_in_complex_estimate=TRUE,
                                 TP=FALSE) {

  TP_result_ids <- as.numeric(unlist(strsplit(manual_stats$TP_result_ids,split=";")))
  FP_result_ids <- as.numeric(unlist(strsplit(manual_stats$FP_result_ids,split=";")))
  FN_manual_ids <- as.numeric(unlist(strsplit(manual_stats$FN_manual_ids,split=";")))

  features <- subset(mapped_features, protein_id == proteinID | protein_id.manual == proteinID)
  #peptides <- unique(unlist(strsplit(features$subunits_annotated, split = ";")))
  peptides <- traces$trace_annotation[protein_id == proteinID]$id
  traces <- subset(traces,trace_ids = peptides)
  traces.long <- toLongFormat(traces$traces)
  setkey(traces.long, id)

  TP_idx=which(features$results_id %in% TP_result_ids)
  FP_idx=which(features$results_id %in% FP_result_ids)
  FN_idx=which(features$manual_id %in% FN_manual_ids)

  features$status = NA
  if (length(TP_idx) > 0){
    features$status[TP_idx]="TP"
  }
  if (length(FP_idx) > 0){
    features$status[FP_idx]="FP"
  }
  if (length(FN_idx) > 0){
    features$status[FN_idx]="FN"
  }

  proteinName = proteinID

  features$status = as.factor(features$status)
  group.colors <- c("TP" = "green3", "FP" = "firebrick1", "FN"="orange")
  group.colors <- group.colors[which(names(group.colors) %in% features$status)]

  grid.newpage()

  if (TP==FALSE){
  p <- ggplot(traces.long) +
        geom_line(aes_string(x='fraction', y='intensity', color='id'), na.rm = TRUE) +
        ggtitle(proteinName) +
        xlab('fraction') +
        ylab('intensity') +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
        theme(plot.margin = unit(c(1,.5,.5,.5),"cm")) +
        theme(plot.title = element_text(vjust=19,size=12)) +
        guides(color=FALSE)
  p <- p + geom_rect(data=features,aes(xmin = left_pp, xmax = right_pp, ymin = 0, ymax = Inf),fill="blue",alpha = 0.15, na.rm = TRUE)
  p <- p + geom_rect(data=features,aes(xmin = left_pp.manual, xmax = right_pp.manual, ymin = 0, ymax = Inf),fill="green3",alpha=0.15, na.rm = TRUE)
  p <- p + geom_vline(data=features,aes(xintercept=apex), linetype="solid",colour="blue",alpha=0.5, na.rm = TRUE)
  p <- p + geom_vline(data=features,aes(xintercept=apex.manual), linetype="solid",colour="green3",alpha=0.5, na.rm = TRUE)
  #print(p)
  }

  if (TP){
  p <- ggplot(traces.long) +
        geom_line(aes_string(x='fraction', y='intensity', color='id'), na.rm = TRUE) +
        ggtitle(proteinName) +
        xlab('fraction') +
        ylab('intensity') +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
        theme(plot.margin = unit(c(1,.5,.5,.5),"cm")) +
        theme(plot.title = element_text(vjust=19, size=12)) +
        guides(color=FALSE)
        p <- p + geom_rect(data=features,aes(xmin = left_pp, xmax = right_pp, ymin = 0, ymax = Inf,fill=status),alpha = 0.15, na.rm = TRUE)
        p <- p + geom_rect(data=features,aes(xmin = left_pp.manual, xmax = right_pp.manual, ymin = 0, ymax = Inf,fill=status),alpha=0.15, na.rm = TRUE)
        p <- p + geom_vline(data=features,aes(xintercept=apex), colour="grey2", linetype="dashed", alpha=0.5,size=1.5, na.rm = TRUE)
        p <- p + geom_vline(data=features,aes(xintercept=apex.manual), colour="green4", linetype="solid",alpha=0.5,size=1.5, na.rm = TRUE)
        p <- p + scale_fill_manual(values=group.colors)
  #print(p)
  }


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
  #dev.off()

}
