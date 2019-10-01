#' Plot pie chart with complex hypothesis completeness
#' @description Plot pie chart showing the number of complex hypotheses with co-elution features of certain completenes.
#' This is only taking into account targets. If decoys are present they are removed using "DECOY" as a flag.
#' @import data.table
#' @param complexFeatures data.table as reported by \code{\link[SECprofiler]{findComplexFeatures}}
#' or \code{\link[SECprofiler]{findProteinFeatures}}.
#' @param hypotheses data.table containing complex hypotheses.
#' Should have the following columns:
#' \itemize{
#' \item complex_id: character strings, a unique id for every complex
#' \item complex_name: character strings, a unique name for every complex
#' \item protein_id: character strings, the protein id, e.g. Uniprot id
#' }
#' @param PDF Logical, whether to generate a PDF file with the summary plot. Default is \code{FALSE}.
#' @param name Character string specifying the name of the PDF file of the pie chart.
#' Only applicable if PDF=\code{TRUE}. Default is "complex_completeness_pie".
#' @examples
#' ## Load example complex feature finding results:
#' complexFeatures <- exampleComplexFeatures
#' complexHypotheses <- exampleComplexHypotheses
#' proteinTraces <- exampleProteinTraces
#' ## Run plotting function:
#' plotSummarizedComplexes(complexFeatures=complexFeatures,
#'                         hypotheses=complexHypotheses,
#'                         protTraces=proteinTraces)
#' @export
plotSummarizedComplexes <- function(complexFeatures,hypotheses,protTraces,PDF=FALSE,name="complex_completeness_pie"){
  targetFeatures <- copy(complexFeatures)
  targetHypotheses <- copy(hypotheses)
  targetFeatures <- targetFeatures[grep("DECOY",targetFeatures$complex_id,invert=TRUE)]
  targetHypotheses <- targetHypotheses[grep("DECOY",targetHypotheses$complex_id,invert=TRUE)]
  proteins_in_targetHypotheses <- unique(targetHypotheses$protein_id)
  proteins_in_traces <- unique(protTraces$traces$id)
  targetHypotheses[,annotated:=1]
  targetHypotheses[,detected:=ifelse(protein_id %in% proteins_in_traces, 1, 0)]
  #targetHypotheses[,protein_collapsed := paste(protein_id,collapse=";"),by=complex_id]
  targetHypotheses[,annotated_collapsed := sum(annotated),by=complex_id]
  targetHypotheses[,detected_collapsed := sum(detected),by=complex_id]
  targetHypotheses[,ms_completeness := detected_collapsed/annotated_collapsed,by=complex_id]
  unique_targetHypotheses <- unique(targetHypotheses,by="complex_id")
  #unique_targetHypotheses_50 <- subset(unique_targetHypotheses,ms_completeness >= 0.5)

  targetFeatures <- getBestFeatures(targetFeatures)
  targetFeatures_min50 <- subset(targetFeatures,(completeness>=0.5) & (completeness<1))
  targetFeatures_100 <- subset(targetFeatures,completeness==1)
  targetFeatures_lower50 <- subset(targetFeatures,completeness<0.5)

  #complexFeatures_noDecoys <- complexFeatures[grep("DECOY",complexFeatures$complex_id,invert=TRUE)]
  complexCompletenessSummary <- data.table(name=c("no co-elution","co-elution\n(< 50% complete)","co-elution\n(>= 50% complete)","co-elution\n(100% complete)"),
                                        count=c(
                                          sum(!(unique_targetHypotheses$complex_id %in% targetFeatures$complex_id)),
                                          sum(unique_targetHypotheses$complex_id %in% targetFeatures_lower50$complex_id),
                                          sum(unique_targetHypotheses$complex_id %in% targetFeatures_min50$complex_id),
                                          sum(unique_targetHypotheses$complex_id %in% targetFeatures_100$complex_id)
                                          )
                                        )


  cbPalette <- c("gray90", "#56B4E9", "#009E73", "royalblue4", "#0072B2", "#D55E00", "#CC79A7")
  if(PDF){pdf(gsub("$|\\.pdf$", ".pdf", name), width = 7, height = 6)}
    print(pie(x=complexCompletenessSummary$count,labels=paste0(complexCompletenessSummary$name,"\n",complexCompletenessSummary$count),col=cbPalette[1:nrow(complexCompletenessSummary)]))
  if(PDF){dev.off()}
}



#' Scatter plot of complex completeness
#' @description Plot the number of observed co-eluting subunits over the number of subunits
#' in the hypothesis, colour-coded by the completeness. This is only taking into account targets.
#' If decoys are present they are removed using "DECOY" as a flag.
#' @import data.table
#' @import ggplot2
#' @param complexFeatures data.table as reported by \code{\link[SECprofiler]{findComplexFeatures}}
#' or \code{\link[SECprofiler]{findProteinFeatures}}.
#' @param PDF Logical, whether to generate a PDF file with the summary plot. Default is \code{FALSE}.
#' @param name Character string specifying the name of the PDF file of the summary plot.
#' Only applicable if PDF=\code{TRUE}. Default is "feature_summary".
#' @examples
#' ## Load example complex feature finding results:
#' complexFeatures <- exampleComplexFeatures
#' ## Run plotting function:
#' plotComplexCompletenessScatter(complexFeatures)
#' @export
plotComplexCompletenessScatter <- function(complexFeatures,PDF=FALSE,name="complex_completeness_scatter"){
  targetFeatures <- copy(complexFeatures)
  targetFeatures <- getBestFeatures(targetFeatures)
  targetFeatures <- targetFeatures[grep("DECOY",targetFeatures$complex_id,invert=TRUE)]
  if(PDF){
    pdf(gsub("$|\\.pdf$", ".pdf", name), width = 5, height = 5)
  }
  p <- ggplot(data=targetFeatures,aes(x=n_subunits_annotated,y=n_subunits_detected,colour=completeness)) +
    geom_point() +
    geom_abline(intercept=0,slope=1) +
    geom_abline(intercept=0,slope=0.5, linetype=2) +
    annotate("text",x=70,y=80, label="100%", angle = 45) +
    annotate("text",x=70,y=9.5, label="50%", angle = 23) +
    scale_x_log10(name = "N subunits in hypothesis", breaks = c(1,10,100), limits = c(1,200)) +
    scale_y_log10(name = "N subunits observed co-eluting", breaks = c(1,10,100), limits = c(1,200)) +
    theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(colour = "Completeness\n") +
    theme(legend.justification=c(0,1), legend.position=c(0.05,0.75)) +
    ggtitle("complex completeness scatter") +
    theme(plot.title = element_text(size=14, face="bold"))
  print(p)
  if(PDF){
    dev.off()
  }
}
