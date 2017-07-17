#' #' Plot the result of the sliding window algorithm.
#' @description Chomatographic visualization
#' @import data.table
#' @import grid
#' @import gridExtra
#' @import gtable
#' @param feature_table data.table
#' @param proteinTraces traces object of type protein
#' @param feature_id character string specifying complex_id
#' @param annotationID character string specifying column name in trace annotation to use for trace labeling, default is protein_id (uniprot id)
#' @param peak_area logical if selected peak area should be highlighted
#' @param sliding_window_are logical if original sliding_window area should be highlighted
#' @param estimated_complex_MW logical if estimated complex MW should be indicated
#' @param monomer_MW logical if monomer MWs should be indicated
#' @param log logical if intensities should be log transformed
#' @export
#' @examples
#' 
#' #------------------------
#' ## Complex level plotting
#' #------------------------
#' 
#' ## Load example data
#' featureTable = exampleComplexFeatures
#' peptideTraces = exampleProteinTraces
#' complexId = "2174-1" # Complex of interest
#' 
#' ## Plot all features found for this Protein
#' plotComplexFeatures(feature_table = featureTable,
#'                     traces = peptideTraces,
#'                     feature_id = complexId,
#'                     peak_area = T,
#'                     onlyBest = F)
#' 
#' #------------------------
#' ## Protein level plotting
#' #------------------------
#'  
#' ## Load example data
#' featureTable = exampleProteinFeatures
#' peptideTraces = examplePeptideTraces
#' proteinId = "P61201" # Protein of interest
#' 
#' ## Plot all features found for this Protein
#' plotComplexFeatures(feature_table = featureTable,
#'                     traces = peptideTraces,
#'                     feature_id = proteinId,
#'                     peak_area = T,
#'                     onlyBest = F)
#' 

plotComplexFeatures <- function(feature_table,
                               traces,
                               feature_id,
                               calibration = NULL,
                               annotationID="protein_id",
                               onlyBest = TRUE,
                               apex=TRUE,
                               peak_area=FALSE,
                               sliding_window_area=FALSE,
                               estimated_complex_MW=FALSE,
                               monomer_MW=FALSE,
                               log=FALSE,
                               legend = TRUE,
                               PDF=FALSE) {
  .tracesTest(traces)
  features <- copy(feature_table)
  if (traces$trace_type == "protein") {
    features <- subset(features, complex_id == feature_id)
    proteins <- unique(unlist(strsplit(features$subunits_with_signal, split = ";")))
    traces <- subset(traces, trace_subset_ids = proteins)
    complexName = unique(features$complex_name)[1]
    n_annotatedSubunits = features$n_subunits_annotated[1]
    n_signalSubunits = features$n_subunits_with_signal[1]
    n_detectedSubunits = features$n_subunits_detected[1]
    complexCompleteness = round(features$completeness[1],digits=2)
    title <- paste0(complexName,
                    "\nAnnotated subunits: ",n_annotatedSubunits,
                    "   Subunits with signal: ",n_signalSubunits,
                    "\nMax. coeluting subunits: ",n_detectedSubunits,
                    "   Max. completeness: ",complexCompleteness)
  } else {
    features <- subset(features, protein_id == feature_id)
    proteins <- unique(unlist(strsplit(features$subunits_annotated, split = ";")))
    traces <- subset(traces, trace_subset_ids = proteins)
    if (annotationID %in% names(traces$trace_annotation)) {
      complexName = traces$trace_annotation[,get(annotationID)][1]
    } else {
      complexName = traces$trace_annotation$protein_id[1]
    }
    n_annotatedSubunits = features$n_subunits_annotated[1]
    n_detectedSubunits = features$n_subunits_detected[1]
    complexCompleteness = round(features$completeness[1],digits=2)
    title <- paste0(complexName,
                    "\nAnnotated subunits: ",n_annotatedSubunits,
                    "\nMax. coeluting subunits: ",n_detectedSubunits,
                    "   Max. completeness: ",complexCompleteness)
  }
  if (onlyBest) {
    features <- getBestFeatures(features)
  }
  traces.long <- toLongFormat(traces$traces)
  traces.long <- merge(traces.long,traces$fraction_annotation,by.x="fraction",by.y="id")
  
  if (annotationID %in% names(traces$trace_annotation)) {
    traces.long <- merge(traces.long,traces$trace_annotation,by="id",all.x=TRUE,all.y=FALSE)
    if (traces$trace_type == "protein") {
      traces.long[,id := get(annotationID)]
      proteins <- traces$trace_annotation[match(proteins,traces$trace_annotation$id)][,get(annotationID)]
    }
    if ("protein_mw" %in% names(traces$trace_annotation)) {
      if (!is.null(calibration)) {
        subunitMW.dt <- data.table(id=traces$trace_annotation[,get(annotationID)],mw=traces$trace_annotation$protein_mw)
        subunitMW.dt$fraction <- calibration$MWtoFraction(subunitMW.dt$mw)
      } else {
        if (monomer_MW){
          message("No calibration. Cannot plot monomer molecular weight.")
          monomer_MW <- FALSE
        }
      }
    } else {
      message("No molecular weight annotation of the traces. Cannot plot monomer molecular weight.")
      monomer_MW <- FALSE
    }
  } else {
    message(paste0(annotationID, " is not present in trace annotation. ID value is used instead."))
    if ("protein_mw" %in% names(traces$trace_annotation)) {
      if (!is.null(calibration)) {
        subunitMW.dt <- data.table(id=traces$trace_annotation$id,mw=traces$trace_annotation$protein_mw)
        subunitMW.dt$fraction <- calibration$MWtoFraction(subunitMW.dt$mw)
      } else {
        if (monomer_MW){
          message("No calibration. Cannot plot monomer molecular weight.")
          monomer_MW <- FALSE
        }
      }
    } else {
      message("No molecular weight annotation of the traces. Cannot plot monomer molecular weight.")
      monomer_MW <- FALSE
    }
  } 
  
  
 
  p <- ggplot(traces.long) +
    geom_line(aes_string(x='fraction', y='intensity', color='id')) +
    xlab('fraction') +
    ylab('intensity') +
    theme_bw() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5),"cm")) +
    ggtitle(title) +
    theme(plot.title = element_text(vjust=19,size=10, face="bold")) +
    guides(fill=FALSE)
  if (log) {
    p <- p + scale_y_log10('log(intensity)')
  }
  if(peak_area==TRUE){
    p <- p + geom_rect(data=features,aes(xmin = left_pp, xmax = right_pp, ymin = -Inf, ymax = Inf), fill="grey", alpha = 0.25)
    p <- p + geom_vline(data=features,aes(xintercept=left_pp), colour="grey",linetype="longdash", alpha=0.4)
    p <- p + geom_vline(data=features,aes(xintercept=right_pp), colour="grey",linetype="longdash", alpha=0.4)
  }
  if(sliding_window_area==TRUE){
    p <- p + geom_vline(data=features,aes(xintercept=left_sw), colour="grey",linetype="longdash", alpha=0.4)
    p <- p + geom_vline(data=features,aes(xintercept=right_sw), colour="grey",linetype="longdash", alpha=0.4)
  }
  if(estimated_complex_MW==TRUE){
    p <- p + geom_vline(data=features,aes(xintercept=complex_sec_estimated), colour="black",linetype="longdash")
  }
  if (monomer_MW==TRUE){
    p <- p + geom_point(data = subunitMW.dt, mapping = aes(x = fraction, y = Inf, colour=id),shape=18,size=5,alpha=.5)
  }
  if (apex==TRUE){
    p <- p + geom_vline(data=features,aes(xintercept=apex), colour="black",linetype="solid")
  }
  if (!legend) {
    p <- p + theme(legend.position="none")
  }
  
  if ("molecular_weight" %in% names(traces$fraction_annotation)) {
    grid.newpage()
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
