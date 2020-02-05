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
      if(!is.null(traces[["genomic_coord"]])){
        traces$genomic_coord <- traces$genomic_coord[trace_ids]
      }
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
#' @param traces Object of class traces.
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
  UseMethod("annotateMolecularWeight", traces)
}

#' @describeIn annotateMolecularWeight Annotate single traces object
#' @export
annotateMolecularWeight.traces <- function(traces, calibration){
  .tracesTest(traces)
  traces$fraction_annotation[,molecular_weight := round(calibration$FractionToMW(traces$fraction_annotation$id),digits=3)]
  traces
}

#' @describeIn annotateMolecularWeight Annotate tracesList object
#' @export
annotateMolecularWeight.tracesList <- function(traces, calibration){
  .tracesListTest(traces)
  res <- lapply(traces, annotateMolecularWeight.traces, calibration = calibration)
  class(res) <- "tracesList"
  .tracesListTest(res)
  return(res)
}

#' Plot traces
#' @description Plot all chromatograms in a traces object. Most generic plotting function.
#' @param traces Object of class traces.
#' @param log Logical, whether the intensities should be plotted in log scale. Default is \code{FALSE}.
#' @param legend Logical, whether a legend of the traces should be plotted. Should be set to \code{FALSE}
#' if many chromatograms are plotted. Default is \code{TRUE}.
#' @param PDF Logical, whether to plot to PDF. PDF file is saved in working directory. Default is \code{FALSE}.
#' @param name Character string with name of the plot, only used if \code{PDF=TRUE}.
#' PDF file is saved under name.pdf. Default is "Traces".
#' @param colorMap named character vector containing valid color specifications for plotting.
#' The names of the vector must correspond to the ids of the peptides to be plotted.
#' @param monomer_MW Logical if monomer MWs should be indicated
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
                        colour_by = "id",
                        highlight=NULL,
                        highlight_col=NULL,
                        colorMap=NULL,
                        monomer_MW=TRUE) {

  .tracesTest(traces)
  traces.long <- toLongFormat(traces$traces)
  traces.long <- merge(traces.long,traces$fraction_annotation,by.x="fraction",by.y="id")
  if(!is.null(highlight)){
    traces.long$outlier <- gsub("\\(.*?\\)","",traces.long$id) %in% gsub("\\(.*?\\)","",highlight)
    if(!any(traces.long$outlier)) highlight <- NULL
  }

  if(colour_by!="id") {
    if(!colour_by %in% names(traces$trace_annotation)){
      stop("colour_by is not availbale in trace_annotation.")
    }
    isoform_annotation <- subset(traces$trace_annotation,select=c("id",colour_by))
    traces.long <- merge(traces.long,isoform_annotation, by.x="id",by.y="id")
    traces.long[,line:=paste0(get(colour_by),id)]
  }

  ## Create a reproducible coloring for the peptides plotted
  if(!is.null(colorMap)){
    if(!all(unique(traces.long$id) %in% names(colorMap))){
      stop("Invalid colorMap specified. Not all traces to be plotted are contained in the colorMap")
    }
  }else{
    cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999")
    ids <- sort(unique(traces.long[[colour_by]]))
    if (length(ids) <= length(cbPalette)) {
      colorMap <- cbPalette[1:length(unique(traces.long[[colour_by]]))]
      names(colorMap) <- ids
    } else {
      colorMap <- createGGplotColMap(unique(traces.long$id))
    }
  }

  if(colour_by == "id") {
    p <- ggplot(traces.long) +
      geom_line(aes_string(x='fraction', y='intensity', colour='id', group='id'))
  } else {

    p <- ggplot(traces.long) +
      geom_line(aes_string(x='fraction', y='intensity', colour=colour_by, group='line'))
  }

  p <- p + xlab('fraction') +
    ylab('intensity') +
    theme_bw() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    theme(plot.margin = unit(c(1,.5,.5,.5),"cm")) +
    ggtitle(name) +
    scale_color_manual(values=colorMap)
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
        geom_line(data = traces.long[outlier == TRUE],
                  aes_string(x='fraction', y='intensity', color='id'), lwd=2) +
        scale_color_manual(values=colorMap, breaks = legend_peps)
      ## scale_color_discrete(breaks = legend_peps)
    }else{
      ## legend_map <- unique(ggplot_build(p)$data[[1]]$colour)
      ## names(legend_map) <- unique(p$data$id)
      ## legend_map[legend_peps] <- highlight_col
      ## legend_vals <- rep(highlight_col, ceiling(length(legend_peps)/ length(highlight_col)))[1:length(legend_peps)]
      p <- p +
        geom_line(data = traces.long[outlier == TRUE],
                  aes_string(x='fraction', y='intensity', lty = 'id'),
                  color = highlight_col, lwd=2)
      # scale_color_discrete(guide = F)
      ## scale_color_manual(values = legend_map, limits = legend_peps)
      # guides(lty = FALSE)
      # scale_color_manual(limits = legend_peps, values = rep(highlight_col, length(legend_peps))) +
      # geom_line(aes_string(x='fraction', y='intensity', color='id'))
    }
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
    if (monomer_MW==TRUE){
      if ("protein_mw" %in% names(traces$trace_annotation)) {
        subunitMW.dt <- data.table(id=traces$trace_annotation$id,mw=traces$trace_annotation$protein_mw)
        subunitMW.dt$fraction <- MWtoFraction(subunitMW.dt$mw)
        subunitMW.dt[,boundary:=MWtoFraction(2*mw)]
        if (length(unique(subunitMW.dt$mw)) > 1) {
          p <- p + geom_point(data = subunitMW.dt, mapping = aes(x = fraction, y = Inf, colour=id),shape=18,size=5,alpha=.5)
        } else {
          p <- p + geom_vline(data = unique(subunitMW.dt), aes(xintercept = fraction), colour="red", linetype="dashed", size=.5)
          p <- p + geom_vline(data = unique(subunitMW.dt), aes(xintercept = boundary), colour="red", linetype="dashed", size=.5, alpha=0.5)
        }
      } else {
        message("No molecular weight annotation of the traces. Cannot plot monomer molecular weight.")
      }
    }
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
  if(plot){
    plot(p)
  }else{
    return(ggplot_gtable(ggplot_build(p)))
  }
  if(PDF){
    dev.off()
  }
}

#' Plot traces
#' @description Plot all chromatograms in a traces object. Most generic plotting function.
#' @param traces Object of class traces.
#' @param log Logical, whether the intensities should be plotted in log scale. Default is \code{FALSE}.
#' @param legend Logical, whether a legend of the traces should be plotted.
#' Should be set to \code{FALSE}
#' if many chromatograms are plotted. Default is \code{TRUE}.
#' @param PDF Logical, whether to plot to PDF. PDF file is saved in working directory.
#'  Default is \code{FALSE}.
#' @param name Character string with name of the plot, only used if \code{PDF=TRUE}.
#' PDF file is saved under name.pdf. Default is "Traces".
#' @param plot Logical, wether to print or return the plot object
#' @param isoformAnnotation Logical, wether to colour traces by their isoform annotation
#' @param colour_by Character string specifying by which column to colour by. Default is id.
#' @param highlight Character vector, ids of the traces to highlight (can be multiple).
#'  Default is \code{NULL}.
#' @param highlight_col Character string, A color to highlight traces in. Must be accepted by ggplot2.
#' @param colorMap named character vector containing valid color specifications for plotting.
#' The names of the vector must correspond to the ids of the peptides to be plotted.
#' @param monomer_MW Logical if monomer MWs should be indicated
#' @export
plot.tracesList <- function(traces,
                            design_matrix = NULL,
                            collapse_conditions = FALSE,
                            aggregateReplicates=FALSE,
                            log=FALSE,
                            legend = TRUE,
                            PDF=FALSE,
                            name="Traces",
                            plot = TRUE,
                            isoformAnnotation = FALSE,
                            colour_by = "id",
                            highlight=NULL,
                            highlight_col=NULL,
                            colorMap=NULL,
                            monomer_MW=TRUE) {
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
  if(colour_by!="id") {
    if(!colour_by %in% names(traces[[1]]$trace_annotation)){
      stop("colour_by is not availbale in trace_annotation.")
    }
    isoform_annotation <- lapply(names(traces), function(tr){subset(traces[[tr]]$trace_annotation,select=c("id",colour_by))})
    isoform_annotation <- unique(do.call("rbind", isoform_annotation))
    traces_long <- merge(traces_long,isoform_annotation, by.x="id",by.y="id")
    traces_long[,line:=paste0(get(colour_by),id)]
  }
  ## Create a common fraction annotation
  traces_frac <- unique(do.call("rbind", lapply(traces, "[[", "fraction_annotation")))
  traces_frac <- unique(subset(traces_frac, select = names(traces_frac) %in% c("id","molecular_weight")))
  traces_long <- merge(traces_long,traces_frac,by.x="fraction",by.y="id")

  if(!is.null(highlight)){
    traces_long$outlier <- gsub("\\(.*?\\)","",traces_long$id) %in% gsub("\\(.*?\\)","",highlight)
    if(!any(traces_long$outlier)) highlight <- NULL
  }

  geom.text.size = 3
  theme.size = (14/5) * geom.text.size

  ## Create a reproducible coloring for the peptides plotted
  if(!is.null(colorMap)){
    if(!all(unique(traces_long[[colour_by]]) %in% names(colorMap))){
      stop("Invalid colorMap specified. Not all traces to be plotted are contained in the colorMap")
    }
  }else{
    cbPalette <- c("#56B4E9","#E69F00", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999")
    ids <- sort(unique(traces_long[[colour_by]]))
    if (length(ids) <= length(cbPalette)) {
      colorMap <- cbPalette[1:length(unique(traces_long[[colour_by]]))]
      names(colorMap) <- ids
    } else {
      colorMap <- createGGplotColMap(ids)
    }
  }

  if (aggregateReplicates){
    traces_long[,meanIntensity := mean(intensity), by=c("Condition","fraction","id")]
    traces_long[,sdIntensity := sd(intensity), by=c("Condition","fraction","id")]
    traces_long <- unique(subset(traces_long,select=c("id","fraction","Condition","molecular_weight","meanIntensity","sdIntensity")))
    traces_long[,Replicate := "average"]
    setnames(traces_long, "meanIntensity", "intensity")
  }

  p <- ggplot(traces_long) +
    xlab('fraction') +
    ylab('intensity') +
    theme_bw() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    theme(plot.margin = unit(c(1,.5,.5,.5),"cm")) +
    ggtitle(name) +
    scale_color_manual(values=colorMap) #+
    #theme(plot.title = element_text(vjust=19,size=10))

  if(collapse_conditions){
    if(colour_by == "id") {
      p <- p + facet_grid(~ Replicate) +
        geom_line(aes_string(x='fraction', y='intensity', color='id', lty = 'Condition'))
      if (aggregateReplicates){
        p <- p + geom_errorbar(aes(x=fraction,ymin=ifelse(intensity-sdIntensity < 0, 0, intensity-sdIntensity), ymax=intensity+sdIntensity, color=id, lty = Condition), width=0.2, size=0.3, position=position_dodge(0.05))
      }
    } else {
      message("Collapsing of conditions is not jet compatible with colouring by your selected id type.
      Plot conditions separately instead.")
      p <- p + facet_grid(Condition ~ Replicate) +
        geom_line(aes_string(x='fraction', y='intensity', color=colour_by, group='line'))
      if (aggregateReplicates){
        p <- p + geom_errorbar(aes(x=fraction,ymin=ifelse(intensity-sdIntensity < 0, 0, intensity-sdIntensity), ymax=intensity+sdIntensity, color=colour_by), width=0.2, size=0.3, position=position_dodge(0.05))
      }
    }
  }else{
    if (length(unique(traces_long$Replicate)) > 1) {
      if(colour_by == "id") {
        p <- p + facet_grid(Condition ~ Replicate) +
          geom_line(aes_string(x='fraction', y='intensity', color='id'))
        if (aggregateReplicates){
          p <- p + geom_errorbar(aes(x=fraction,ymin=ifelse(intensity-sdIntensity < 0, 0, intensity-sdIntensity), ymax=intensity+sdIntensity, color=id), width=0.2, size=0.3, position=position_dodge(0.05))
        }
      } else {
        p <- p + facet_grid(Condition ~ Replicate) +
          geom_line(aes_string(x='fraction', y='intensity', color=colour_by, group='line'))
      }
      if (aggregateReplicates){
        p <- p + geom_errorbar(aes(x=fraction,ymin=ifelse(intensity-sdIntensity < 0, 0, intensity-sdIntensity), ymax=intensity+sdIntensity, color=colour_by), width=0.2, size=0.3, position=position_dodge(0.05))
      }
    } else {
      if(colour_by == "id") {
        p <- p + facet_grid(Condition ~ .) +
          geom_line(aes_string(x='fraction', y='intensity', color='id'))
        if (aggregateReplicates){
          p <- p + geom_errorbar(aes(x=fraction,ymin=ifelse(intensity-sdIntensity < 0, 0, intensity-sdIntensity), ymax=intensity+sdIntensity, color=id), width=0.2, size=0.3, position=position_dodge(0.05))
        }
      } else {
        p <- p + facet_grid(Condition ~ .) +
          geom_line(aes_string(x='fraction', y='intensity', color=colour_by, group='line'))
        if (aggregateReplicates){
          p <- p + geom_errorbar(aes(x=fraction,ymin=ifelse(intensity-sdIntensity < 0, 0, intensity-sdIntensity), ymax=intensity+sdIntensity, color=colour_by), width=0.2, size=0.3, position=position_dodge(0.05))
        }
      }
    }
  }

  if(!is.null(highlight)){
    legend_peps <- unique(traces_long[outlier == TRUE, id])
    if(is.null(highlight_col)){
      if(collapse_conditions){
        p <- p +
          geom_line(data = traces_long[outlier == TRUE], aes_string(x='fraction', y='intensity', color='id', lty = 'Condition'), lwd=2) +
          scale_color_manual(values = colorMap, breaks = legend_peps)
      }else{
        p <- p +
          geom_line(data = traces_long[outlier == TRUE], aes_string(x='fraction', y='intensity', color='id'), lwd=2) +
          scale_color_manual(values = colorMap, breaks = legend_peps)
      }

    }else{
      ## legend_map <- unique(ggplot_build(p)$data[[1]]$colour)
      ## names(legend_map) <- unique(p$data$id)
      ## legend_map[legend_peps] <- highlight_col
      ## legend_vals <- rep(highlight_col, ceiling(length(legend_peps)/ length(highlight_col)))[1:length(legend_peps)]
      if(collapse_conditions){
        p <- p +
          geom_line(data = traces_long[outlier == TRUE],
                    aes(x=fraction, y=intensity, lty = Condition, group = interaction(Condition, id), color = id),
                    lwd=2) +
          # scale_color_discrete(guide = F)
          scale_color_manual(values = colorMap, breaks = legend_peps)
        # guides(lty = FALSE)
        # scale_color_manual(limits = legend_peps, values = rep(highlight_col, length(legend_peps))) +
        # geom_line(aes_string(x='fraction', y='intensity', color='id'))
      }else{
        p <- p +
          geom_line(data = traces_long[outlier == TRUE], aes_string(x='fraction', y='intensity', color = 'id'),
                    lwd=2) +
          # scale_color_discrete(guide = F)
          scale_color_manual(values = colorMap, breaks = legend_peps)
        ## scale_color_manual(values = legend_map, limits = legend_peps)
        # guides(lty = FALSE)
        # scale_color_manual(limits = legend_peps, values = rep(highlight_col, length(legend_peps))) +
        # geom_line(aes_string(x='fraction', y='intensity', color='id'))
      }

    }
  }

  if ("molecular_weight" %in% names(traces_frac)) {
    fraction_ann <- traces_frac
    tr <- lm(log(fraction_ann$molecular_weight) ~ fraction_ann$id)
    intercept <- as.numeric(tr$coefficients[1])
    slope <- as.numeric(tr$coefficients[2])
    mwtransform <- function(x){exp(slope*x + intercept)}
    MWtoFraction <- function(x){round((log(x)-intercept)/(slope), digits = 0)}
    mw <- round(fraction_ann$molecular_weight, digits = 0)
    breaks_MW <- mw[seq(1,length(mw), length.out = length(seq(min(traces_frac$id),
                                                              max(traces_frac$id),10)))]
    p <- p + scale_x_continuous(name="fraction",
                                breaks=seq(min(traces_frac$id),
                                           max(traces_frac$id),10),
                                labels=seq(min(traces_frac$id),
                                           max(traces_frac$id),10),
                                sec.axis = dup_axis(trans = ~.,
                                                    breaks=seq(min(traces_frac$id),
                                                               max(traces_frac$id),10),
                                                    labels = breaks_MW,
                                                    name = "MW (kDa)"))
    if (monomer_MW==TRUE){
      if ("protein_mw" %in% names(traces[[1]]$trace_annotation)) {
        ann_tab <- lapply(traces, function(t){subset(t$trace_annotation, select=c(colour_by,"protein_mw"))})
        ann_tab <- unique(do.call(rbind,ann_tab))
        subunitMW.dt <- data.table(id=ann_tab[[colour_by]],mw=ann_tab$protein_mw)
        subunitMW.dt$fraction <- MWtoFraction(subunitMW.dt$mw)
        subunitMW.dt[,boundary:=MWtoFraction(2*mw)]
        if (length(unique(subunitMW.dt$mw)) > 1) {
          p <- p + geom_point(data = subunitMW.dt, mapping = aes(x = fraction, y = Inf, colour=id), shape=18,size=5,alpha=.5)
        } else {
          p <- p + geom_vline(data = unique(subunitMW.dt), aes(xintercept = fraction), colour="red", linetype="dashed", size=.5)
          p <- p + geom_vline(data = unique(subunitMW.dt), aes(xintercept = boundary), colour="red", linetype="dashed", size=.5, alpha=0.5)
        }
      } else {
        message("No molecular weight annotation of the traces. Cannot plot monomer molecular weight.")
      }
    }
  } else {
    p <- p + scale_x_continuous(name="fraction",
                                breaks=seq(min(traces$fraction_annotation$id),
                                           max(traces$fraction_annotation$id),10),
                                labels=seq(min(traces$fraction_annotation$id),
                                           max(traces$fraction_annotation$id),10))
  }

  if (log) {
    p <- p + scale_y_log10('log(intensity)')
  }

  p <- p + theme(axis.text = element_text(size = theme.size, colour="black"))

  p <- p + theme(legend.position="bottom") +
    theme(legend.text=element_text(size=theme.size)) +
    theme(legend.title=element_blank()) +
    guides(col = guide_legend(ncol = 4))

  if (!legend) {
    p <- p + theme(legend.position="none")
  }

  if (length(ids) > 40) {
    p <- p + theme(legend.position="none")
  }

  if(PDF){
    pdf(paste0(name,".pdf"),width=5,height=4)
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

  if(traces$trace_type == "peptide"){
    no_ptraces <- length(unique(traces$trace_annotation$protein_id))
    no_pdecoys <- length(unique(grep("DECOY", traces$trace_annotation$protein_id)))
    no_ptargets <- no_ptraces - no_pdecoys
    pct_pdecoys <- signif(no_pdecoys/no_ptraces * 100, 2)
    pres <- c(no_ptraces, no_ptargets, no_pdecoys, pct_pdecoys)
    names(pres) <- c("No. of Proteins", "No. of Targets", "No. of Decoys", "% Decoys")
    res <- list(Peptides = res, Proteins = pres)
  }
  annotation_info <- names(traces$trace_annotation)
  fraction_info = length(traces$fraction_annotation$id)
  type=traces$trace_type

  summary=list(metrics=res,type=type,annotations=annotation_info,fraction_count=fraction_info)
  if("SibPepCorr" %in% names(traces$trace_annotation)) {
    SibPepCorr_summary <- summary(traces$trace_annotation$SibPepCorr)
    summary$SibPepCorr_summary <- SibPepCorr_summary
  }
  if("RepPepCorr" %in% names(traces$trace_annotation)) {
    RepPepCorr_summary <- summary(traces$trace_annotation$RepPepCorr)
    summary$RepPepCorr_summary <- RepPepCorr_summary
  }
  summary
}

#' Summarize a tracesList object
#' @description Summarize a tracesList object to get an overview.
#' @param traces Object of class tracesList
#' @return A summary list with following entries:
#'          \itemize{
#'           \item \code{metrics} Summary of traces counts including decoy statistics.
#'           \item \code{type} Type of traces: peptides or proteins.
#'           \item \code{annotations} Names of trace annotations.
#'           \item \code{fraction_count} Number of fractions in traces object.
#'           \item \code{SibPepCorr} Summary of sibling peptide correlation.
#'           Only reported if sibling peptide correlations were previously calculated.
#'          }
#' @export

summary.tracesList <- function(traces) {
  .tracesListTest(traces)
  res <- lapply(names(traces),function(tr){
    trace <- traces[[tr]]
    cat(paste0("###########################\n## ",
           tr, "\n",
           "###########################\n\n"))
    print(summary.traces(trace))
  })

}

#' Summarize a tracesList object
#' @description Summarize a tracesList object to get an overview.
#' @param traces Object of class tracesList
#' @return A summary list with following entries:
#'          \itemize{
#'           \item \code{metrics} Summary of traces counts including decoy statistics.
#'           \item \code{type} Type of traces: peptides or proteins.
#'           \item \code{annotations} Names of trace annotations.
#'           \item \code{fraction_count} Number of fractions in traces object.
#'           \item \code{SibPepCorr} Summary of sibling peptide correlation.
#'           Only reported if sibling peptide correlations were previously calculated.
#'          }
#' @export
print.tracesList <- function(traces){
  summary.tracesList(traces)
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

print.traces <- function(traces){
  summary.traces(traces)
}

#' Test if an object is of class traces.
#' @param traces Object of class traces.
#' @param type Character string specifying whether a specific type of traces is required.
#' @param additionalItems Character string specifying additional entries that are required in the list.
#' The two options are "peptide" or "protein". Default is \code{NULL},
#' meaning that no specific type is required.
.tracesTest <- function(traces,type=NULL, additionalItems=NULL){
  if (! class(traces)=="traces") {
    stop("Object is not of class traces.")
  }
  if (! all(c("traces","trace_type","trace_annotation","fraction_annotation") %in% names(traces))) {
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
  if(!is.null(additionalItems)){
    contained <- (additionalItems %in% names(traces))
    if(!all(contained)){
      stop(paste0("Required entries not found: ", additionalItems[!contained]))
    }
  }
}

#' Test if an object is of class tracesList.
#' @param traces Object of class tracesList.
#' @param type Character string specifying whether a specific type of traces is required.
#' @param additionalItems Character string specifying additional entries that are required in the list.
#' The two options are "peptide" or "protein". Default is \code{NULL},
#' meaning that no specific type is required.
.tracesListTest <- function(tracesList, type=NULL, additionalItems=NULL){
  if (! class(tracesList)=="tracesList") {
    stop("Object is not of class tracesList")
  }
  if(is.null(names(tracesList))) stop("TracesList must consist of named traces objects. No names detected.")
  res <- lapply(tracesList, function(traces){
    if (! all(c("traces","trace_type","trace_annotation","fraction_annotation") %in% names(traces))) {
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
    if(!is.null(additionalItems)){
      contained <- (additionalItems %in% names(traces))
      if(!all(contained)){
        stop(paste0("Required entries not found: ", additionalItems[!contained]))
      }
    }
  })
}

#' Update trace and fraction annotation
#' @description  Add information to trace and fraction annotation
#' @param traces Object of class traces.
#' @return Object of class traces.
#' @export
#'
updateTraces <- function(traces) {
  UseMethod("updateTraces", traces)
}

#' @describeIn updateTraces Update trace and fraction annotation
#' @export
updateTraces.traces <- function(traces) {
  # call all functions that compute summary statistics that could have changed due to the manipulation
  # of the traces object (note that these functions should be relatively fast)
  traces <- annotateFractions(traces)
  .tracesTest(traces)
  return(traces)
}

#' @describeIn updateTraces Update trace and fraction annotation
#' @export
updateTraces.tracesList <- function(traces) {
  traces <- annotateFractions(traces)
  .tracesListTest(traces)
  return(traces)
}
