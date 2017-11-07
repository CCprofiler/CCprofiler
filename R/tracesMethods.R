# Due to: http://stackoverflow.com/questions/24501245/data-table-throws-object-not-found-error
.datatable.aware=TRUE

#' Subset traces or fractions in traces object
#' @description  Subset a taces object by a specified column in trace_annotation and/or fraction ids.
#' @param traces Object of class traces.
#' @param trace_subset_ids Character vector specifying the trace identifiers
#'        for subsetting traces, e.g. peptide or protein ids.
#' @param trace_subset_type Character string  specifying the column name
#'        for applying the trace_subset_ids filter, defaults to "id".
#' @param fraction_ids Numeric vector specifying the fraction identifiers
#'        for subsetting traces. The resulting traces object will not have
#'        fraction ids starting from 1, co-elution feature finding might thus
#'        be compromised. Please don't use this option for general
#'        pre-processing purposes.
#' @return Object of class traces.
#' @examples
#' # Load some example data:
#' inputPeptideTraces <- examplePeptideTraces
#' inputProteinTraces <- exampleProteinTraces
#' subsetPeptides <- c("AIIDEFEQK","AIQLSGAEQLEALK","AKEALIAASETLK")
#' subsetProtein <- "Q15021"
#' subsetFractions <- c(5:20)
#' # Run subsetting and inspect resulting traces object:
#' peptideSubsettedPeptideTraces <- subset(inputPeptideTraces,trace_subset_ids=subsetPeptides,fraction_ids=subsetFractions)
#' summary(peptideSubsettedPeptideTraces)
#'
#' proteinSubsettedPeptideTraces <- subset(inputPeptideTraces,trace_subset_ids=subsetProtein,trace_subset_type="protein_id")
#' summary(proteinSubsettedPeptideTraces)
#'
#' proteinSubsettedProteinTraces <- subset(inputProteinTraces,trace_subset_ids=subsetProtein)
#' summary(proteinSubsettedProteinTraces)
#' @export
#' 

subset.traces <- function(traces,trace_subset_ids=NULL,trace_subset_type="id",fraction_ids=NULL){
  .tracesTest(traces)
  if (!is.null(trace_subset_ids)) {
    if (trace_subset_type %in% names(traces$trace_annotation)) {
      traces$trace_annotation <- subset(traces$trace_annotation, get(trace_subset_type) %in% trace_subset_ids)
      trace_ids <- traces$trace_annotation$id
      traces$traces <- subset(traces$traces,id %in% trace_ids)
      if (nrow(traces$traces) == 0) {
        message("Caution! Subsetting returns empty traces object.")
      }
    } else {
      stop(paste0(trace_subset_type, "is not a valid trace_subset_type."))
    }
  }
  if (!is.null(fraction_ids)){
    if(class(fraction_ids) == "numeric"){
      fraction_ids <- as.character(fraction_ids)
    }
    fraction_ids <- intersect(names(traces$traces), fraction_ids)
    traces$traces <- subset(traces$traces,select=c(fraction_ids,"id"))
    traces$fraction_annotation <- subset(traces$fraction_annotation ,id %in% fraction_ids)
    if (nrow(traces$traces) == 0) {
      stop(paste0("fraction_ids (",fraction_ids,") do not match the available fractions in the traces object."))
    }
  }
  traces
}

#' Subset traces or fractions in tracesList object
#' @description  Subset a taces object by a specified column in trace_annotation and/or fraction ids.
#' @param traces Object of class traces.
#' @param trace_subset_ids Character vector specifying the trace identifiers
#'        for subsetting traces, e.g. peptide or protein ids.
#' @param trace_subset_type Character string  specifying the column name
#'        for applying the trace_subset_ids filter, defaults to "id".
#' @param fraction_ids Numeric vector specifying the fraction identifiers
#'        for subsetting traces. The resulting traces object will not have
#'        fraction ids starting from 1, co-elution feature finding might thus
#'        be compromised. Please don't use this option for general
#'        pre-processing purposes.
#' @return Object of class traces.
#' @export
#' 

subset.tracesList <- function(tracesList,
                              trace_subset_ids=NULL,
                              trace_subset_type="id",
                              fraction_ids=NULL){
  .tracesListTest(tracesList)
  res <- lapply(tracesList, subset.traces,
                trace_subset_ids=trace_subset_ids,
                trace_subset_type=trace_subset_type,
                fraction_ids=fraction_ids)
  class(res) <- "tracesList"
  .tracesListTest(res)
  return(res)
}

#' Get intensity matrix from traces object
#' @description  Get a matrix of intensity values from a traces object.
#' @param traces Object of class traces.
#' @return Numeric matrix with intensity values.
#' @examples
#' intensityMatrix <- getIntensityMatrix(examplePeptideTraces)
#' head(intensityMatrix)
#' @export
getIntensityMatrix <- function(traces){
  .tracesTest(traces)
  tr <- copy(traces$traces)
  ids <- tr$id
  intensity.mat <- data.matrix(tr[,id := NULL])
  rownames(intensity.mat) <- ids
  intensity.mat
}


#' Convert traces from wide format to long format
#' @description Convert the data.table traces from a trace object from wide format to long format.
#' @param traces.dt A data.table with an id column \code{id} and
#'        columns of continuously numbered fractions.
#' @return A data.table with columns
#'          \itemize{
#'           \item \code{id}
#'           \item \code{fraction}
#'           \item \code{intensity}
#'          }
#' @examples
#' tracesWide <- examplePeptideTraces$traces
#' tracesLong <- toLongFormat(tracesWide)
#' tracesLong
#' @export
toLongFormat <- function(traces.dt) {
  traces.dt.long <-
    melt(traces.dt, id.var='id', variable.name='fraction',
         value.name='intensity', variable.factor=FALSE)
  traces.dt.long[, fraction := as.numeric(fraction)]
  setkey(traces.dt.long,id)
  data.table(traces.dt.long)
  traces.dt.long
}

#' Annotate traces with molecular weight calibration.
#' @description Annotate fractions in traces object with calibrated molecular weight.
#' This only applies to specific fractionation strategies e.g. SEC)
#' @param Object of class traces.
#' @return Object of class traces with moleclar weight annotation of fractions.
#' @examples
#' # Load relevant data:
#' inputTraces <- examplePeptideTraces
#' calibrationTable <- exampleCalibrationTable
#' # Perform molecular weight calibration:
#' calibration = calibrateMW(calibrationTable)
#' # Perform molecular weight annotation:
#' mwTraces <- annotateMolecularWeight(inputTraces, calibration)
#' @export
annotateMolecularWeight <- function(traces, calibration){
  .tracesTest(traces)
  traces$fraction_annotation[,molecular_weight := round(calibration$FractionToMW(traces$fraction_annotation$id),digits=3)]
  traces
}

#' Plot traces
#' @description Plot all chromatograms in a traces object. Most generic plotting function.
#' @param Object of class traces.
#' @param log Logical, whether the intensities should be plotted in log scale. Default is \code{FALSE}.
#' @param legend Logical, whether a legend of the traces should be plotted. Should be set to \code{FALSE}
#' if many chromatograms are plotted. Default is \code{TRUE}.
#' @param PDF Logical, whether to plot to PDF. PDF file is saved in working directory. Default is \code{FALSE}.
#' @param name Character string with name of the plot, only used if \code{PDF=TRUE}.
#' PDF file is saved under name.pdf. Default is "Traces".
#' @examples
#' # Protein traces
#' proteinTraces=exampleProteinTraces
#' plot(proteinTraces)
#' # Annotate traces with molecular weight to include molecular weight information in plot.
#' calibrationTable <- exampleCalibrationTable
#' # Perform molecular weight calibration:
#' calibration = calibrateMW(calibrationTable)
#' # Perform molecular weight annotation:
#' mwProteinTraces <- annotateMolecularWeight(proteinTraces, calibration)
#' plot(mwProteinTraces)
#' # Plot all peptides of a specific protein across a defined chromatographic region
#' peptideTraces <- examplePeptideTraces
#' subsetPeptideTraces <- subset(peptideTraces,trace_subset_ids="Q9UHV9",trace_subset_type="protein_id",fraction_ids=c(30:70))
#' plot(subsetPeptideTraces,legend=FALSE)
#' @export
plot.traces <- function(traces,
                        log=FALSE,
                        legend = TRUE,
                        PDF=FALSE,
                        name="Traces",
                        plot = TRUE,
                        highlight=NULL,
                        highlight_col=NULL) {
  
  .tracesTest(traces)
  traces.long <- toLongFormat(traces$traces)
  traces.long <- merge(traces.long,traces$fraction_annotation,by.x="fraction",by.y="id")
  if(!is.null(highlight)){
    traces.long$outlier <- gsub("\\(.*?\\)","",traces.long$id) %in% gsub("\\(.*?\\)","",highlight)
    if(!any(traces.long$outlier)) highlight <- NULL
  }
  
  p <- ggplot(traces.long) +
    geom_line(aes_string(x='fraction', y='intensity', color='id')) +
    xlab('fraction') +
    ylab('intensity') +
    theme_bw() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    theme(plot.margin = unit(c(1,.5,.5,.5),"cm")) +
    ggtitle(name)#+
  #theme(plot.title = element_text(vjust=19,size=10))
  if (log) {
    p <- p + scale_y_log10('log(intensity)')
  }
  if (!legend) {
    p <- p + theme(legend.position="none")
  }
  if(!is.null(highlight)){
    legend_peps <- unique(traces.long[outlier == TRUE, id])
    if(is.null(highlight_col)){
      p <- p + 
        geom_line(data = traces.long[outlier == TRUE], aes_string(x='fraction', y='intensity', color='id'), lwd=2) +
        scale_color_discrete(breaks = legend_peps)
    }else{
      legend_map <- unique(ggplot_build(p)$data[[1]]$colour)
      names(legend_map) <- unique(p$data$id)
      legend_map[legend_peps] <- highlight_col
      legend_vals <- rep(highlight_col, ceiling(length(legend_peps)/ length(highlight_col)))[1:length(legend_peps)]
      p <- p + 
        geom_line(data = traces.long[outlier == TRUE], aes_string(x='fraction', y='intensity', lty = 'id'), color = highlight_col, lwd=2) +
        # scale_color_discrete(guide = F)
        scale_color_manual(values = legend_map, limits = legend_peps) 
      # guides(lty = FALSE)
      # scale_color_manual(limits = legend_peps, values = rep(highlight_col, length(legend_peps))) +
      # geom_line(aes_string(x='fraction', y='intensity', color='id'))
    }
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
    if(plot){
      grid.draw(g)
    }else{
      return(g)
    }
    
    if(PDF){
      dev.off()
    }
  }else{
    if(PDF){
      pdf(paste0(name,".pdf"))
    }
    if(plot){
      plot(p)
    }else{
      return(ggplot_gtable(ggplot_build(p)))
    }
    
    if(PDF){
      dev.off()
    }
  }
}

#' Plot traces
#' @description Plot all chromatograms in a traces object. Most generic plotting function.
#' @param Object of class traces.
#' @param log Logical, whether the intensities should be plotted in log scale. Default is \code{FALSE}.
#' @param legend Logical, whether a legend of the traces should be plotted. 
#' Should be set to \code{FALSE}
#' if many chromatograms are plotted. Default is \code{TRUE}.
#' @param PDF Logical, whether to plot to PDF. PDF file is saved in working directory.
#'  Default is \code{FALSE}.
#' @param name Character string with name of the plot, only used if \code{PDF=TRUE}.
#' PDF file is saved under name.pdf. Default is "Traces".
#' @param plot Logical, wether to print or return the plot object
#' @param highlight Character vector, ids of the traces to highlight (can be multiple).
#'  Default is \code{NULL}.
#' @param highlight_col Character string, A color to highlight traces in. Must be accepted by ggplot2.
#' @export

plot.tracesList <- function(traces,
                            design_matrix = NULL,
                            collapse_conditions = FALSE,
                            log=FALSE,
                            legend = TRUE,
                            PDF=FALSE,
                            name="Traces",
                            plot = TRUE,
                            highlight=NULL,
                            highlight_col=NULL) {
  .tracesListTest(traces)
  if(!is.null(design_matrix)){
    if(!all(design_matrix$Sample_name %in% names(traces))){
      stop("Invalid design matrix")
    }
  }else{
    design_matrix <- data.table(Sample_name = names(traces),
                                Condition = "", 
                                Replicate = 1:length(traces))
  }
  tracesList <- lapply(names(traces), function(tr){
    res <- toLongFormat(traces[[tr]]$traces)
    res$Condition <- design_matrix[Sample_name == tr, Condition]
    res$Replicate <- design_matrix[Sample_name == tr, Replicate]
    res
  })
  traces_long <- do.call("rbind", tracesList)
  traces_frac <- unique(do.call("rbind", lapply(traces, "[[", "fraction_annotation")))
  traces_long <- merge(traces_long,traces_frac,by.x="fraction",by.y="id")
  if(!is.null(highlight)){
    traces_long$outlier <- gsub("\\(.*?\\)","",traces_long$id) %in% gsub("\\(.*?\\)","",highlight)
    if(!any(traces_long$outlier)) highlight <- NULL
  }
  
  p <- ggplot(traces_long) +
    xlab('fraction') +
    ylab('intensity') +
    theme_bw() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    theme(plot.margin = unit(c(1,.5,.5,.5),"cm")) +
    ggtitle(name) +
    theme(plot.title = element_text(vjust=19,size=10))
  if(collapse_conditions){
    p <- p + facet_grid(~ Replicate) +
      geom_line(aes_string(x='fraction', y='intensity', color='id', lty = 'Condition'))
  }else{
    p <- p + facet_grid(Condition ~ Replicate) +
      geom_line(aes_string(x='fraction', y='intensity', color='id'))
    
  }
  if(!is.null(highlight)){
    legend_peps <- unique(traces_long[outlier == TRUE, id])
    if(is.null(highlight_col)){
      if(collapse_conditions){
        p <- p + 
          geom_line(data = traces_long[outlier == TRUE], aes_string(x='fraction', y='intensity', color='id', lty = 'Condition'), lwd=2) +
          scale_color_discrete(breaks = legend_peps)
      }else{
        p <- p + 
          geom_line(data = traces_long[outlier == TRUE], aes_string(x='fraction', y='intensity', color='id'), lwd=2) +
          scale_color_discrete(breaks = legend_peps)
      }
      
    }else{
      legend_map <- unique(ggplot_build(p)$data[[1]]$colour)
      names(legend_map) <- unique(p$data$id)
      legend_map[legend_peps] <- highlight_col
      legend_vals <- rep(highlight_col, ceiling(length(legend_peps)/ length(highlight_col)))[1:length(legend_peps)]
      if(collapse_conditions){
        p <- p + 
          geom_line(data = traces_long[outlier == TRUE], 
                    aes(x=fraction, y=intensity, lty = Condition, group = interaction(Condition, id), color = id),
                     lwd=2) +
          # scale_color_discrete(guide = F)
          scale_color_manual(values = legend_map, limits = legend_peps) 
        # guides(lty = FALSE)
        # scale_color_manual(limits = legend_peps, values = rep(highlight_col, length(legend_peps))) +
        # geom_line(aes_string(x='fraction', y='intensity', color='id'))
      }else{
        p <- p + 
          geom_line(data = traces_long[outlier == TRUE], aes_string(x='fraction', y='intensity', lty = 'id', color = 'id'),
                    lwd=2) +
          # scale_color_discrete(guide = F)
          scale_color_manual(values = legend_map, limits = legend_peps) 
        # guides(lty = FALSE)
        # scale_color_manual(limits = legend_peps, values = rep(highlight_col, length(legend_peps))) +
        # geom_line(aes_string(x='fraction', y='intensity', color='id'))
      }
      
    }
  }
  
  if (log) {
    p <- p + scale_y_log10('log(intensity)')
  }
  if (!legend) {
    p <- p + theme(legend.position="none")
  }
  if(PDF){
    pdf(paste0(name,".pdf"))
  }
  if(plot){
    plot(p)
  }else{
    return(p)
  }
  if(PDF){
    dev.off()
  }
}

#' Summarize a traces object
#' @description Summarize a traces object to get an overview.
#' @param traces Object of class traces.
#' @return A summary list with following entries:
#'          \itemize{
#'           \item \code{metrics} Summary of traces counts including decoy statistics.
#'           \item \code{type} Type of traces: peptides or proteins.
#'           \item \code{annotations} Names of trace annotations.
#'           \item \code{fraction_count} Number of fractions in traces object.
#'           \item \code{SibPepCorr} Summary of sibling peptide correlation.
#'           Only reported if sibling peptide correlations were previously calculated.
#'          }
#' @examples
#' summary(examplePeptideTracesFiltered)
#' @export
summary.traces <- function(traces) {
  .tracesTest(traces)
  no_traces <- nrow(traces$trace_annotation)
  no_decoys <- length(grep("DECOY", traces$trace_annotation$protein_id))
  no_targets <- no_traces - no_decoys
  pct_decoys <- signif(no_decoys/no_traces * 100, 2)
  res <- c(no_traces, no_targets, no_decoys, pct_decoys)
  names(res) <- c("No. of Traces", "No. of Targets", "No. of Decoys", "% Decoys")
  annotation_info <- names(traces$trace_annotation)
  fraction_info = length(traces$fraction_annotation$id)
  type=traces$trace_type
  if("SibPepCorr" %in% names(traces$trace_annotation)) {
    SibPepCorr_summary <- summary(traces$trace_annotation$SibPepCorr)
    summary=list(metrics=res,type=type,annotations=annotation_info,fraction_count=fraction_info,SibPepCorr_summary=SibPepCorr_summary)
  } else {
    summary=list(metrics=res,type=type,annotations=annotation_info,fraction_count=fraction_info)
  }
  summary
}

#' Test if an object is of class traces.
#' @param traces Object of class traces.
#' @param type Character string specifying whether a specific type of traces is required.
#' The two options are "peptide" or "protein". Default is \code{NULL},
#' meaning that no specific type is required.
.tracesTest <- function(traces,type=NULL){
  if (! class(traces)=="traces") {
    stop("Object is not of class traces.")
  }
  if (! all(names(traces)==c("traces","trace_type","trace_annotation","fraction_annotation"))) {
    stop("Traces object doesn't contain all necessary items: traces, trace_type, trace_annotation, and fraction_annotation.")
  }
  if (!is.null(type)) {
    if (type != traces$trace_type) {
      stop("Traces object is of wrong type. Please check your input traces.")
    }
  }
  if (! identical(traces$traces$id,traces$trace_annotation$id)) {
    stop("IDs in traces and trace_annotation are not identical.")
  }
  if (! identical(names(traces$traces),c(traces$fraction_annotation$id,"id"))) {
    stop("Fractions in traces and fraction_annotation are not identical.")
  }
}

#' Test if an object is of class tracesList.
#' @param traces Object of class tracesList.
#' @param type Character string specifying whether a specific type of traces is required.
#' The two options are "peptide" or "protein". Default is \code{NULL},
#' meaning that no specific type is required.
.tracesListTest <- function(tracesList,type=NULL){
  if (! class(tracesList)=="tracesList") {
    stop("Object is not of class tracesList")
  }
  if(is.null(names(tracesList))) stop("TracesList must consist of named traces objects. No names detected.")
  res <- lapply(tracesList, function(traces){
    if (! all(names(traces)==c("traces","trace_type","trace_annotation","fraction_annotation"))) {
      stop("At least one traces object doesn't contain all necessary items: traces, trace_type, trace_annotation, and fraction_annotation.")
    }
    if (!is.null(type)) {
      if (type != traces$trace_type) {
        stop("At least one traces object is of wrong type. Please check your input traces.")
      }
    }
    if (! identical(traces$traces$id,traces$trace_annotation$id)) {
      stop("In at least one traces object: IDs in traces and trace_annotation are not identical.")
    }
    if (! identical(names(traces$traces),c(traces$fraction_annotation$id,"id"))) {
      stop("In at least one traces object: Fractions in traces and fraction_annotation are not identical.")
    }
    
  })
}

