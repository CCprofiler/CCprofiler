#' Plot MS coverage of the hypotheses in the PCP data
#' @description Generate plots that describe the coverage of the complex hypotheses in the PCP data.
#' @import data.table
#' @import ggplot2
#' @param hypotheses data.table containing complex hypotheses.
#' Should have the following columns:
#' \itemize{
#' \item complex_id: character strings, a unique id for every complex
#' \item complex_name: character strings, a unique name for every complex
#' \item protein_id: character strings, the protein id, e.g. Uniprot id
#' }
#' @param protTraces An object of class traces, type="protein".
#' @param PDF Logical whether to save plots in PDF, default = \code{FALSE}.
#' @return Multiple plots summarizing the MS coverage.
#' @export
plotSummarizedMScoverage <- function(hypotheses,protTraces,PDF=FALSE, name_suffix="complexHypotheses"){
  .tracesTest(protTraces, type="protein")
  hypotheses <- hypotheses[grep("DECOY",hypotheses$complex_id,invert=TRUE)]
  proteins_in_hypotheses <- unique(hypotheses$protein_id)
  proteins_in_traces <- unique(protTraces$traces$id)
  proteinMScoverageSummary <- data.table(name=c("not detected","detected"),count=c(sum(!(proteins_in_hypotheses %in% proteins_in_traces)),sum(proteins_in_hypotheses %in% proteins_in_traces)))
  plotProteinMScoveragePie(proteinMScoverageSummary,PDF=PDF,name_suffix=name_suffix)
  hypotheses[,annotated:=1]
  hypotheses[,detected:=ifelse(protein_id %in% proteins_in_traces, 1, 0)]
  #hypotheses[,protein_collapsed := paste(protein_id,collapse=";"),by=complex_id]
  hypotheses[,annotated_collapsed := sum(annotated),by=complex_id]
  hypotheses[,detected_collapsed := sum(detected),by=complex_id]
  hypotheses[,ms_completeness := detected_collapsed/annotated_collapsed,by=complex_id]
  unique_hypotheses <- unique(hypotheses,by="complex_id")
  complexMScoverageSummary <- data.table(name=c("not detected","< 50% detected",">= 50% detected","fully detected"),
                                         count=c(sum(unique_hypotheses$ms_completeness == 0),
                                                 sum((unique_hypotheses$ms_completeness > 0) & (unique_hypotheses$ms_completeness < 0.5)),
                                                 sum((unique_hypotheses$ms_completeness < 1) & (unique_hypotheses$ms_completeness >= 0.5)),
                                                 sum(unique_hypotheses$ms_completeness == 1))
  )
  plotComplexMScoveragePie(complexMScoverageSummary,PDF=PDF,name_suffix=name_suffix)
  plotComplexMScoverageScatter(unique_hypotheses,PDF=PDF,name_suffix=name_suffix)
}

#' plotProteinMScoveragePie.
#' @description plotProteinMScoveragePie
#' @import data.table
#' @import ggplot2
#' @param proteinMScoverageSummary data.table with summary
#' @param PDF logical default = TRUE
#' @return plots summarizing the MS coverage
plotProteinMScoveragePie <- function(proteinMScoverageSummary,PDF=FALSE, name_suffix="complexHypotheses"){
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  if(PDF){
    pdf(paste0("proteinMScoveragePie_",name_suffix,".pdf"))
  }
  print(pie(x=proteinMScoverageSummary$count,
            labels=paste0(proteinMScoverageSummary$name,"\n",proteinMScoverageSummary$count),
            col=cbPalette[1:nrow(proteinMScoverageSummary)],
            main = paste0("Protein coverage ",name_suffix)))
  if(PDF){
    dev.off()
  }
}

#' plotComplexMScoveragePie.
#' @description plotComplexMScoveragePie
#' @import data.table
#' @import ggplot2
#' @param complexMScoverageSummary data.table with summary
#' @param PDF logical default = TRUE
#' @return plots summarizing the MS coverage
plotComplexMScoveragePie <- function(complexMScoverageSummary,PDF=TRUE, name_suffix="complexHypotheses"){
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  if(PDF){
    pdf(paste0("complexMScoveragePie_",name_suffix,".pdf"))
  }
  print(pie(x=complexMScoverageSummary$count,
            labels=paste0(complexMScoverageSummary$name,"\n",complexMScoverageSummary$count),
            col=cbPalette[1:nrow(complexMScoverageSummary)],
            main = paste0("Hypothesis coverage ", name_suffix)))
  if(PDF){
    dev.off()
  }
}

#' plotComplexMScoverageScatter.
#' @description plotComplexMScoverageScatter
#' @import data.table
#' @import ggplot2
#' @param complexMScoverage data.table with summary
#' @param PDF logical default = TRUE
#' @return plots summarizing the MS coverage
plotComplexMScoverageScatter <- function(complexMScoverage,PDF=FALSE, name_suffix="complexHypotheses"){
  if(PDF){
    pdf(paste0("complexMScoverageScatter_",name_suffix,".pdf"))
  }
  p <- ggplot(data=complexMScoverage,aes(x=annotated_collapsed,y=detected_collapsed,colour=ms_completeness)) +
    geom_point() +
    geom_abline(intercept=0,slope=1) +
    geom_abline(intercept=0,slope=0.5, linetype=2) +
    annotate("text",x=70,y=80, label="100%", angle = 45) +
    annotate("text",x=70,y=9.5, label="50%", angle = 23) +
    scale_x_log10(name = "N subunits in hypothesis", breaks = c(1,10,100), limits = c(1,200)) +
    scale_y_log10(name = "N subunits observed by MS", breaks = c(1,10,100), limits = c(1,200)) +
    theme_classic() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(colour = "MS completeness\n") +
    theme(legend.justification=c(0,1), legend.position=c(0.05,0.75)) +
    ggtitle(paste0("Hypothesis coverage ",name_suffix)) + 
    theme(plot.title = element_text(size=14, face="bold"))
  print(p)
  if(PDF){
    dev.off()
  }
}
