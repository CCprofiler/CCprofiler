#' Plot protein features.
#'
#' @param res proteinFeatures
#' @param pepTraces An object of type \code{traces.obj}.
#' @param proteinID protein ID
#' @param window.size Size of the window. Numeric.
#' @param plot_peak logical
#' @param plot_monomer logical
#' @export

plot.proteinFeatures <- function(res,
                                 traces.obj,
                                 proteinID,
                                 plot_peak=TRUE,
                                 plot_monomer=TRUE,
                                 log=FALSE,
                                 plot_apex=TRUE,
                                 plot_in_complex_estimate=TRUE) {


  features <- subset(res, protein_id == proteinID)
  #peptides <- unique(unlist(strsplit(features$subunits_annotated, split = ";")))
  peptides <- traces.obj$trace_annotation[protein_id == proteinID]$id
  traces <- subset(traces.obj,trace_ids = peptides)
  traces.long <- toLongFormat(traces$traces)

  # add monomer MW
  monomer = as.numeric(unlist(features$monomer_sec))
  monomer_max = as.numeric(unlist(features$monomer_sec_max))

  setkey(traces.long, id)

  proteinName = unique(features$protein_id)[1]

  group.colors <- c("TRUE" = "green3", "FALSE" = "firebrick1")

  p <- ggplot(traces.long) +
    geom_line(aes_string(x='fraction', y='intensity', color='id')) +
    ggtitle(proteinName) +
    xlab('fraction') +
    ylab('intensity') +
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
    p <- p + geom_vline(data=features,aes(xintercept=apex), colour="darkgrey", linetype="solid")
  }
  #p <- p + theme(legend.position="none")
  print(p)
}
