#' Plot pie chart with complex hypothesis completeness
#' @description Plot pie chart showing the number of complex hypotheses with co-elution features of certain completenes.
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
#' @export
plotSummarizedComplexes <- function(complexFeatures,hypotheses,protTraces,PDF=FALSE,name="complex_completeness_pie"){
  proteins_in_hypotheses <- unique(hypotheses$protein_id)
  proteins_in_traces <- unique(protTraces$traces$id)
  hypotheses[,annotated:=1]
  hypotheses[,detected:=ifelse(protein_id %in% proteins_in_traces, 1, 0)]
  #hypotheses[,protein_collapsed := paste(protein_id,collapse=";"),by=complex_id]
  hypotheses[,annotated_collapsed := sum(annotated),by=complex_id]
  hypotheses[,detected_collapsed := sum(detected),by=complex_id]
  hypotheses[,ms_completeness := detected_collapsed/annotated_collapsed,by=complex_id]
  unique_hypotheses <- unique(hypotheses,by="complex_id")
  unique_hypotheses_50 <- subset(unique_hypotheses,ms_completeness >= 0.5)

  complexFeatures <- getBestComplexFeature(complexFeatures)
  complexFeatures_min50 <- subset(complexFeatures,completeness>=0.5)
  complexFeatures_lower50 <- subset(complexFeatures,completeness<0.5)


  complexFeatures_noDecoys <- complexFeatures[grep("DECOY",complexFeatures$complex_id,invert=TRUE)]
  complexCompletenessSummary <- data.table(name=c("no co-elution","co-elution\n(< 50% complete)","co-elution\n(>= 50% complete)"),
                                        count=c(
                                          sum(!(unique_hypotheses_50$complex_id %in% complexFeatures$complex_id)),
                                          sum(unique_hypotheses_50$complex_id %in% complexFeatures_lower50$complex_id),
                                          sum(unique_hypotheses_50$complex_id %in% complexFeatures_min50$complex_id)
                                          )
                                        )


  cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  if(PDF){pdf(gsub("$|\\.pdf$", ".pdf", name))}
    pie(x=complexCompletenessSummary$count,labels=paste0(complexCompletenessSummary$name,"\n",complexCompletenessSummary$count),col=cbPalette[1:nrow(complexCompletenessSummary)])
  if(PDF){dev.off()}
}



#' Scatter plot of complex completeness
#' @description Plot the number of observed co-eluting subunits over the number of subunits 
#' in the hypothesis, colour-coded by the completeness.
#' @import data.table
#' @import ggplot2
#' @param complexFeatures data.table as reported by \code{\link[SECprofiler]{findComplexFeatures}} 
#' or \code{\link[SECprofiler]{findProteinFeatures}}.
#' @param PDF Logical, whether to generate a PDF file with the summary plot. Default is \code{FALSE}.
#' @param name Character string specifying the name of the PDF file of the summary plot.
#' Only applicable if PDF=\code{TRUE}. Default is "feature_summary".
#' @export
plotComplexCompletenessScatter <- function(complexFeatures,PDF=FALSE,name="complex_completeness_scatter"){
  complexFeatures <- getBestComplexFeature(complexFeatures)
  complexFeatures <- complexFeatures[grep("DECOY",complexFeatures$complex_id,invert=TRUE)]
  if(PDF){
    pdf(gsub("$|\\.pdf$", ".pdf", name))
  }
  p <- ggplot(data=complexFeatures,aes(x=n_subunits_annotated,y=n_subunits_detected,colour=completeness)) +
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


