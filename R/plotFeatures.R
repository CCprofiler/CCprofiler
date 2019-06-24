#' #' #' Plot the result of the sliding window algorithm.
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
                                colour_by = "id",
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
                                colour_by = "id",
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
      if(is.character(features$completeness)){
        complexCompleteness = features$completeness
      } else {
        complexCompleteness = round(features$completeness[1],digits=2)
      }
      title <- paste0(complexName,
                      "\nAnnotated subunits: ",n_annotatedSubunits,
                      "   Subunits with signal: ",n_signalSubunits,
                      "\nMax. coeluting subunits: ",n_detectedSubunits,
                      "   Max. completeness: ",complexCompleteness)
    } else {
      title=paste0(features$complex_name,"\n ","\n ","\n ")
    }
  } else {
    if ("protein_id" %in% names(features)) {
      features <- subset(features, protein_id == feature_id)
    } else if ("proteoform_id" %in% names(features)) {
      features <- subset(features, proteoform_id == feature_id)
    } else {
      stop("Features are not protein or proteoform level.")
    }
    proteins <- unique(unlist(strsplit(features$subunits_annotated, split = ";")))
    traces <- subset(traces, trace_subset_ids = proteins)
    if (annotation_label %in% names(traces$trace_annotation)) {
      complexName = traces$trace_annotation[,get(annotation_label)][1]
    } else {
      if ("protein_id" %in% names(features)) {
        complexName = traces$trace_annotation$protein_id[1]
      } else if ("proteoform_id" %in% names(features)) {
        complexName = traces$trace_annotation$proteoform_id[1]
      }
    }
    n_annotatedSubunits = features$n_subunits_annotated[1]
    n_detectedSubunits = features$n_subunits_detected[1]
    if(is.character(features$completeness)){
      complexCompleteness = features$completeness[1]
    } else {
      complexCompleteness = round(features$completeness[1],digits=2)
    }
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
  if(colour_by!="id") {
    if(!colour_by %in% names(traces$trace_annotation)){
      stop("colour_by is not availbale in trace_annotation.")
    }
    isoform_annotation <- subset(traces$trace_annotation,select=c("id",colour_by))
    traces.long <- merge(traces.long,isoform_annotation, by.x="id",by.y="id")
    traces.long[,line:=paste0(get(colour_by),id)]
  }

  ## Get traces to highlight
  if(!is.null(highlight)){
    traces.long$outlier <- gsub("\\(.*?\\)","",traces.long$id) %in% gsub("\\(.*?\\)","",highlight)
  }

  if (annotation_label %in% names(traces$trace_annotation)) {
    traces.long <- merge(traces.long,traces$trace_annotation,by=intersect(names(traces.long),names(traces$trace_annotation)),all.x=TRUE,all.y=FALSE)
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
    if(!all(unique(traces_long$id) %in% names(colorMap))){
      stop("Invalid colorMap specified. Not all traces to be plotted are contained in the colorMap")
    }
  }else{
    cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999")
    ids <- sort(unique(traces.long[[colour_by]]))
    if (length(ids) <= length(cbPalette)) {
      colorMap <- cbPalette[1:length(unique(traces.long[[colour_by]]))]
      names(colorMap) <- ids
    } else {
      colorMap <- createGGplotColMap(ids)
    }
  }

  p <- ggplot(traces.long) +
    #geom_line(aes_string(x='fraction', y='intensity', color='id')) +
    xlab('fraction') +
    ylab('intensity') +
    theme_bw() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5),"cm")) +
    ggtitle(title) +
    theme(plot.title = element_text(vjust=19,size=10, face="bold")) +
    guides(fill=FALSE) +
    scale_color_manual(values=colorMap)

  if(colour_by == "id") {
    p <- p + geom_line(aes_string(x='fraction', y='intensity', color='id'))
  } else {
    p <- p + geom_line(aes_string(x='fraction', y='intensity', color=colour_by, group='line'))
  }


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
    if (length(unique(subunitMW.dt$id)) > 1) {
      p <- p + geom_point(data = subunitMW.dt, mapping = aes(x = fraction, y = Inf, colour=id),shape=18,size=5,alpha=.5)
    } else {
      p <- p + geom_vline(data = unique(subunitMW.dt), aes(xintercept = fraction), colour="red", linetype="dashed", size=1)
    }
  }
  if (apex==TRUE){
    p <- p + geom_vline(data=features,aes(xintercept=apex), colour="black",linetype="solid")
  }
  if (!legend) {
    p <- p + theme(legend.position="none")
  }

  if ("molecular_weight" %in% names(traces$fraction_annotation)) {
    fraction_ann <- traces$fraction_annotation
    tr <- lm(log(fraction_ann$molecular_weight) ~ fraction_ann$id)
    intercept <- as.numeric(tr$coefficients[1])
    slope <- as.numeric(tr$coefficients[2])
    mwtransform <- function(x){exp(slope*x + intercept)}
    MWtoFraction <- function(x){round((log(x)-intercept)/(slope), digits = 0)}
    mw <- round(fraction_ann$molecular_weight, digits = 0)
    breaks_MW <- mw[seq(1,length(mw), length.out = length(seq(min(traces$fraction_annotation$id),
                                                              max(traces$fraction_annotation$id),10)))]
    p <- p + scale_x_continuous(name="fraction",
                                breaks=seq(min(traces$fraction_annotation$id),
                                           max(traces$fraction_annotation$id),10),
                                labels=seq(min(traces$fraction_annotation$id),
                                           max(traces$fraction_annotation$id),10),
                                sec.axis = dup_axis(trans = ~.,
                                                    breaks=seq(min(traces$fraction_annotation$id),
                                                               max(traces$fraction_annotation$id),10),
                                                    labels = breaks_MW,
                                                    name = "MW (kDa)"))
  } else {
    p <- p + scale_x_continuous(name="fraction",
                                breaks=seq(min(traces$fraction_annotation$id),
                                           max(traces$fraction_annotation$id),10),
                                labels=seq(min(traces$fraction_annotation$id),
                                           max(traces$fraction_annotation$id),10))
  }

  if(PDF){
    pdf(paste0(name,".pdf"))
  }
  plot(p)
  if(PDF){
    dev.off()
  }

}
