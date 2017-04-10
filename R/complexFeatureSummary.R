#' plotSummarizedComplexes
#' @description plotSummarizedComplexes
#' @import data.table
#' @import ggplot2
#' @param complexFeatures data.table with detected complexes
#' @param hypotheses data.table with complex hypotheses
#' @param PDF logical default = TRUE
#' @return plots summarizing the complex subfeatures
#' @export
plotSummarizedComplexes <- function(complexFeatures,hypotheses,protTraces,PDF=TRUE){
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
  if(PDF){pdf("complexCompletenessPie.pdf")}
    print(pie(x=complexCompletenessSummary$count,labels=paste0(complexCompletenessSummary$name,"\n",complexCompletenessSummary$count),col=cbPalette[1:nrow(complexCompletenessSummary)]))
  if(PDF){dev.off()}
}

#' plotComplexSubfeatureSummary
#' @description plotComplexSubfeatureSummary
#' @import data.table
#' @import ggplot2
#' @param complexFeatures data.table with detected complexes
#' @param PDF logical default = TRUE
#' @return plots summarizing the complex subfeatures
#' @export
plotComplexSubfeatureSummary <- function(complexFeatures,PDF=TRUE){
  complexFeatures <- complexFeatures[grep("DECOY",complexFeatures$complex_id,invert=TRUE)]
  complexFeatures[,n_features := sum(.N), by="complex_id"]
  complexFeatures <- unique(complexFeatures, by="complex_id")

  if(PDF){
    pdf("complexFeatureSummary.pdf")
  }
  p <- ggplot(data=complexFeatures,aes(x=n_features)) +
    stat_bin(binwidth=1) +
    stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-1.5) +
    labs(x="number of features",y="number of complexes") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
  print(p)
  if(PDF){
    dev.off()
  }
}

#' plotComplexCompletenessScatter
#' @description plotComplexCompletenessScatter
#' @import data.table
#' @import ggplot2
#' @param complexFeatures data.table with detected complexes
#' @param PDF logical default = TRUE
#' @return plots summarizing the completeness
#' @export
plotComplexCompletenessScatter <- function(complexFeatures,PDF=TRUE){
  complexFeatures <- getBestComplexFeature(complexFeatures)
  complexFeatures <- complexFeatures[grep("DECOY",complexFeatures$complex_id,invert=TRUE)]
  if(PDF){
    pdf("complexCompletenessScatter.pdf")
  }
  p <- ggplot(data=complexFeatures,aes(x=n_subunits_annotated,y=n_subunits_detected,colour=completeness)) +
    geom_point() +
    geom_abline(intercept=0,slope=1) +
    geom_abline(intercept=0,slope=0.5, linetype=2) +
    annotate("text",x=70,y=80, label="100%", angle = 45) +
    annotate("text",x=70,y=9.5, label="50%", angle = 23) +
    scale_x_log10(name = "N subunits in hypothesis", breaks = c(1,10,100), limits = c(1,200)) +
    scale_y_log10(name = "N subunits observed co-eluting", breaks = c(1,10,100), limits = c(1,200)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
    labs(colour = "Completeness\n") +
    theme(legend.justification=c(0,1), legend.position=c(0.05,0.75))
  print(p)
  if(PDF){
    dev.off()
  }
}

#' Summarize complex features
#' @description Summarize complex features
#' @export
summarizeComplexFeatures <- function(res){
  total_confirmed_hypotheses = length(unique(res$complex_id))
  total_features = nrow(res)
  sub=res
  sub[ , `:=`( COUNT = .N , IDX = 1:.N ) , by = complex_id ]
  sub = subset(sub,COUNT>=2)
  setkey(sub, "complex_id")
  total_hypotheses_with_multiple_features = nrow(unique(sub))
  data.table(total_confirmed_hypotheses = total_confirmed_hypotheses,
    total_features = total_features,
    total_hypotheses_with_multiple_features = total_hypotheses_with_multiple_features)
}

#' Visual complex feature summary
#' @description Visual complex feature summary
#' @export
plotComplexFeatureSummary <- function(res){
  subunit_density = ggplot(res) +
    geom_density(aes(x=n_subunits_detected))
  subunit_histogram = ggplot(res) +
    geom_histogram(aes(x=n_subunits_detected),binwidth = 1)
  completeness_density = ggplot(res) +
    geom_density(aes(x=completeness))
  completeness_histogram = ggplot(res) +
    geom_histogram(aes(x=completeness),binwidth = 0.05)
  peak_corr_density = ggplot(res) +
    geom_density(aes(x=peak_corr))
  peak_corr_histogram = ggplot(res) +
    geom_histogram(aes(x=peak_corr),binwidth = 0.05)
  multiplot(subunit_density, completeness_density, peak_corr_density, subunit_histogram, completeness_histogram, peak_corr_histogram, cols=2)
}

## http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
