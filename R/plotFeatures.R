#' #' Plot the result of the sliding window algorithm.
#' @description Chomatographic visualization
#' @import data.table
#' @import grid
#' @import gridExtra
#' @import gtable
#' @param feature_table data.table
#' @param proteinTraces traces object of type protein
#' @param feature_id Character string specifying complex_id
#' @param annotation_label Character string specifying column name in trace annotation to use for trace labeling, default is protein_id (uniprot id).
#' @param peak_area Logical if selected peak area should be highlighted.
#' @param sliding_window_are Logical if original sliding_window area should be highlighted
#' @param estimated_complex_MW Logical if estimated complex MW should be indicated
#' @param monomer_MW Logical if monomer MWs should be indicated
#' @param log Logical if intensities should be log transformed, default is \code{fFALSE}.
#' @param colorMap named character vector containing valid color specifications for plotting.
#' The names of the vector must correspond to the ids of the peptides to be plotted.
#' @export
#' @examples
#'
#' #------------------------
#' ## Complex level plotting
#' #------------------------
#'
#' ## Load example data
#' featureTable = exampleComplexFeatures
#' proteinTraces = exampleProteinTraces
#' complexId = "2174-1" # Complex of interest
#'
#' ## Plot all features found for this Protein
#' plotFeatures(feature_table = featureTable,
#'                     traces = proteinTraces,
#'                     feature_id = complexId,
#'                     peak_area = TRUE,
#'                     onlyBest = FALSE)
#'
#' #------------------------
#' ## Protein level plotting
#' #------------------------
#'
#' ## Load example data
#'  featureTable = exampleProteinFeatures
#'  peptideTraces = examplePeptideTraces
#'  proteinId = "P61201" # Protein of interest
#'
#'  ## Plot all features found for this Protein
#'  plotFeatures(feature_table = featureTable,
#'                      traces = peptideTraces,
#'                      feature_id = proteinId,
#'                      peak_area = TRUE,
#'                      onlyBest = FALSE,
#'                      legend = FALSE)
#'

plotFeatures <- function(feature_table,
                                traces,
                                feature_id,
                                calibration = NULL,
                                annotation_label="protein_id",
                                onlyBest = TRUE,
                                apex=TRUE,
                                peak_area=FALSE,
                                sliding_window_area=FALSE,
                                estimated_complex_MW=FALSE,
                                highlight=NULL,
                                highlight_col=NULL,
                                monomer_MW=FALSE,
                                log=FALSE,
                                legend = TRUE,
                                PDF=FALSE,
                                name = "Traces",
                                colorMap=NULL) {
  UseMethod("plotFeatures", traces)
}

plotFeatures.traces <- function(feature_table,
                         traces,
                         feature_id,
                         calibration = NULL,
                         annotation_label="protein_id",
                         onlyBest = TRUE,
                         apex=TRUE,
                         peak_area=FALSE,
                         sliding_window_area=FALSE,
                         estimated_complex_MW=FALSE,
                         highlight=NULL,
                         highlight_col=NULL,
                         monomer_MW=FALSE,
                         log=FALSE,
                         legend = TRUE,
                         PDF=FALSE,
                         name = "Traces",
                         colorMap=NULL) {
  .tracesTest(traces)
  features <- copy(feature_table)
  if (traces$trace_type == "protein") {
    features <- subset(features, complex_id == feature_id)
    if ("subunits_with_signal" %in% names(features)){
      proteins <- unique(unlist(strsplit(features$subunits_with_signal, split = ";")))
    } else {
      proteins <- unique(unlist(strsplit(features$subunits, split = ";")))
    }
    traces <- subset(traces, trace_subset_ids = proteins)
    complexName = unique(features$complex_name)[1]
    if ("complexCompleteness" %in% names(features)) {
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
      title=paste0(features$complex_name,"\n ","\n ","\n ")
    }
  } else {
    features <- subset(features, protein_id == feature_id)
    proteins <- unique(unlist(strsplit(features$subunits_annotated, split = ";")))
    traces <- subset(traces, trace_subset_ids = proteins)
    if (annotation_label %in% names(traces$trace_annotation)) {
      complexName = traces$trace_annotation[,get(annotation_label)][1]
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

  ## Get traces to highlight
  if(!is.null(highlight)){
    traces.long$outlier <- gsub("\\(.*?\\)","",traces.long$id) %in% gsub("\\(.*?\\)","",highlight)
  }

  if (annotation_label %in% names(traces$trace_annotation)) {
    traces.long <- merge(traces.long,traces$trace_annotation,by="id",all.x=TRUE,all.y=FALSE)
    if (traces$trace_type == "protein") {
      traces.long[,id := get(annotation_label)]
      proteins <- traces$trace_annotation[match(proteins,traces$trace_annotation$id)][,get(annotation_label)]
    }
    if ("protein_mw" %in% names(traces$trace_annotation)) {
      if (!is.null(calibration)) {
        subunitMW.dt <- data.table(id=traces$trace_annotation[,get(annotation_label)],mw=traces$trace_annotation$protein_mw)
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
    message(paste0(annotation_label, " is not present in trace annotation. ID value is used instead."))
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

  ## Create a reproducible coloring for the peptides plotted
  if(!is.null(colorMap)){
    if(!all(unique(traces.long$id) %in% names(colorMap))){
      stop("Invalid colorMap specified. Not all traces to be plotted are contained in the colorMap")
    }
  }else{
    colorMap <- createGGplotColMap(unique(traces.long$id))
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
    guides(fill=FALSE) +
    scale_color_manual(values=colorMap)

  if(!is.null(highlight)){
    legend_peps <- unique(traces.long[outlier == TRUE, id])
    if(is.null(highlight_col)){
      p <- p +
        geom_line(data = traces.long[outlier == TRUE], aes_string(x='fraction', y='intensity', color='id'), lwd=2) +
        scale_color_manual(values=colorMap, breaks = legend_peps)
    }else{
      p <- p +
        geom_line(data = traces.long[outlier == TRUE], aes_string(x='fraction', y='intensity', group = 'id'), color = highlight_col, lwd=2)
      # scale_color_discrete(breaks = legend_peps)
    }
  }

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

#' plotFeatures.tracesList
#' @param feature_table data.table
#' @param traces traces object of type peptide or protein.
#' @param feature_id Character string specifying complex_id
#' @param design_matrix data.table with design matrix
#' @param calibration calibration function
#' @param annotation_label Character string specifying column name in trace annotation to use for trace labeling, default is protein_id (uniprot id).
#' @param peak_area Logical if selected peak area should be highlighted.
#' @param onlyBest logical true or false
#' @param apex logical true or false
#' @param peak_area logical true or false
#' @param sliding_window_area Logical if original sliding_window area should be highlighted
#' @param estimated_complex_MW Logical if estimated complex MW should be indicated
#' @param monomer_MW Logical if monomer MWs should be indicated
#' @param log Logical if intensities should be log transformed, default is \code{fFALSE}.
#' @param colorMap named character vector containing valid color specifications for plotting.
#' The names of the vector must correspond to the ids of the peptides to be plotted.
#' @export
plotFeatures.tracesList <- function(feature_table,
                               traces,
                               feature_id,
                               design_matrix = NULL,
                               calibration = NULL,
                               annotation_label="protein_id",
                               onlyBest = TRUE,
                               apex=TRUE,
                               peak_area=FALSE,
                               sliding_window_area=FALSE,
                               estimated_complex_MW=FALSE,
                               highlight=NULL,
                               highlight_col=NULL,
                               monomer_MW=FALSE,
                               log=FALSE,
                               legend = TRUE,
                               PDF=FALSE,
                               name = "Traces",
                               colorMap=NULL) {
  ## .tracesListTest(traces)
  features <- copy(feature_table)

  if(!is.null(design_matrix)){
    if(!all(design_matrix$Sample_name %in% names(traces))){
      stop("Invalid design matrix")
    }
  }else{
    design_matrix <- data.table(Sample_name = names(traces),
                                Condition = "",
                                Replicate = 1:length(traces))
  }


  if (traces[[1]]$trace_type == "protein") {
    features <- subset(features, complex_id == feature_id)
    if ("subunits_with_signal" %in% names(features)){
      proteins <- unique(unlist(strsplit(features$subunits_with_signal, split = ";")))
    } else {
      proteins <- unique(unlist(strsplit(features$subunits, split = ";")))
    }
    traces <- subset(traces, trace_subset_ids = proteins)
    traceAnn <- do.call(rbind, lapply(traces, "[[", "trace_annotation"))
    traceAnn <- unique(traceAnn,by="id")
    complexName = unique(features$complex_name)[1]
    if ("complexCompleteness" %in% names(features)) {
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
      title=paste0(features$complex_name,"\n ","\n ","\n ")
    }
  } else {
    features <- subset(features, protein_id == feature_id)
    proteins <- unique(unlist(strsplit(features$subunits_annotated, split = ";")))
    traces <- subset(traces, trace_subset_ids = proteins)
    traceAnn <- do.call(rbind, lapply(traces, "[[", "trace_annotation"))
    traceAnn <- unique(traceAnn)
    if (annotation_label %in% names(traceAnn)) {
      complexName = traceAnn[,get(annotation_label)][1]
    } else {
      complexName = traceAnn$protein_id[1]
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
  tracesList <- lapply(names(traces), function(tr){
    res <- toLongFormat(traces[[tr]]$traces)
    res$Condition <- design_matrix[Sample_name == tr, Condition]
    res$Replicate <- design_matrix[Sample_name == tr, Replicate]
    res
  })
  traces_long <- do.call("rbind", tracesList)
  ## traces.long <- toLongFormat(traces$traces)
  ## traces.long <- merge(traces.long,traces$fraction_annotation,by.x="fraction",by.y="id")

  ## Get traces to highlight
  if(!is.null(highlight)){
    traces_long$outlier <- gsub("\\(.*?\\)","",traces_long$id) %in% gsub("\\(.*?\\)","",highlight)
  }

  ## Generate the plot labels
  if (annotation_label %in% names(traceAnn)) {
    traces_long <- merge(traces_long, traceAnn, by="id",all.x=TRUE,all.y=FALSE)
    if (traces[[1]]$trace_type == "protein") {
      traces_long[,id := get(annotation_label)]
      proteins <- traceAnn[match(proteins,traceAnn$id)][,get(annotation_label)]
    }
    if ("protein_mw" %in% names(traceAnn)) {
      if (!is.null(calibration)) {
        subunitMW.dt <- data.table(id=traceAnn[,get(annotation_label)],mw=traceAnn$protein_mw)
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
    message(paste0(annotation_label, " is not present in trace annotation. ID value is used instead."))
    if ("protein_mw" %in% names(traceAnn)) {
      if (!is.null(calibration)) {
        subunitMW.dt <- data.table(id=traceAnn[,get(annotation_label)],mw=traceAnn$protein_mw)
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

  ## Create a reproducible coloring for the peptides plotted
  if(!is.null(colorMap)){
    if(!all(unique(traces_long$id) %in% names(colorMap))){
      stop("Invalid colorMap specified. Not all traces to be plotted are contained in the colorMap")
    }
  }else{
    colorMap <- createGGplotColMap(unique(traces_long$id))
  }

  p <- ggplot(traces_long) +
    geom_line(aes_string(x='fraction', y='intensity', color='id')) +
    xlab('fraction') +
    ylab('intensity') +
    theme_bw() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5),"cm")) +
    ggtitle(title) +
    theme(plot.title = element_text(vjust=15,size=10, face="bold")) +
    scale_color_manual(values=colorMap) +
    guides(fill=FALSE)

  ## Split the plot with design matrix
  if (length(unique(traces_long$Replicate)) > 1) {
    p <- p + facet_grid(Condition ~ Replicate)
  } else {
    p <- p + facet_grid(Condition ~ .)
  }
  ## Highlight traces
  if(!is.null(highlight)){
    legend_peps <- unique(traces_long[outlier == TRUE, id])
    if(is.null(highlight_col)){
      p <- p +
        geom_line(data = traces_long[outlier == TRUE], aes_string(x='fraction', y='intensity', color='id'), lwd=2) +
        scale_color_manual(values=colorMap, breaks = legend_peps)
        ## scale_color_discrete(breaks = legend_peps)
    }else{
      p <- p +
        geom_line(data = traces_long[outlier == TRUE], aes_string(x='fraction', y='intensity', group = 'id'), color = highlight_col, lwd=2)
      # scale_color_discrete(breaks = legend_peps)
    }
  }

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

  if ("molecular_weight" %in% names(traces[[1]]$fraction_annotation)) {
    message("Molecular weight annotation for multiple conditions not yet implemented")
    ## grid.newpage()
    ## p2 <- p
    ## p <- p + scale_x_continuous(name="fraction",
    ##                             breaks=seq(min(traces[[1]]$fraction_annotation$id),
    ##                                        max(traces[[1]]$fraction_annotation$id),10),
    ##                             labels=seq(min(traces[[1]]$fraction_annotation$id),
    ##                                        max(traces[[1]]$fraction_annotation$id),10))
    ## p2 <- p2 + scale_x_continuous(name="molecular weight (kDa)",
    ##                               breaks=seq(min(traces[[1]]$fraction_annotation$id),max(traces[[1]]$fraction_annotation$id),10),
    ##                               labels=round(traces[[1]]$fraction_annotation$molecular_weight,digits=0)[seq(1,length(traces[[1]]$fraction_annotation$id),10)]
    ## )
    ## ## extract gtable
    ## g1 <- ggplot_gtable(ggplot_build(p))
    ## g2 <- ggplot_gtable(ggplot_build(p2))
    ## ## overlap the panel of the 2nd plot on that of the 1st plot
    ## pp <- c(subset(g1$layout, name=="panel", se=t:r))

    ## g <- gtable_add_grob(g1
    ## g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name=="panel")]], pp$t, pp$l, pp$b, pp$l)
    ## ## steal axis from second plot and modify
    ## ia <- which(g2$layout$name == "axis-b")
    ## ga <- g2$grobs[[ia]]
    ## ax <- ga$children[[2]]
    ## ## switch position of ticks and labels
    ## ax$heights <- rev(ax$heights)
    ## ax$grobs <- rev(ax$grobs)
    ## ## modify existing row to be tall enough for axis
    ## g$heights[[2]] <- g$heights[g2$layout[ia,]$t]
    ## ## add new axis
    ## g <- gtable_add_grob(g, ax, 2, 4, 2, 4)
    ## ## add new row for upper axis label
    ## g <- gtable_add_rows(g, g2$heights[1], 1)
    ## ## steal axis label from second plot
    ## ia2 <- which(g2$layout$name == "xlab-b")
    ## ga2 <- g2$grobs[[ia2]]
    ## g <- gtable_add_grob(g,ga2, 3, 4, 2, 4)

    ## if(PDF){
    ##   pdf(paste0(name,".pdf"))
    ## }
    ## grid.draw(g)
    ## if(PDF){
    ##   dev.off()
    ## }
  ## }else{
  }
  if(PDF){
    pdf(paste0(name,".pdf"))
  }
  plot(p)
  if(PDF){
    dev.off()
  }
}


#' Create a reproducible color map for a set of ids by ordering
#' them alphabetically
#' @importFrom scales hue_pal
createGGplotColMap <- function(ids){
  ids <- sort(unique(ids))
  colormap <- hue_pal()(length(ids))
  names(colormap) <- ids
  colormap
}
