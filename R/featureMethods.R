#' Summarize co-elution features
#' @description Summarize co-elution features
#' @param feature_table data.table as reported by \code{\link[SECprofiler]{findComplexFeatures}} 
#' or \code{\link[SECprofiler]{findProteinFeatures}}.
#' @param plot Logical, whether to generate a feature summary plot. Default is \code{TRUE}.
#' @param PDF Logical, whether to generate a PDF file with the summary plot. Default is \code{FALSE}.
#' @param name Character string specifying the name of the PDF file of the summary plot.
#' Only applicable if PDF=\code{TRUE}. Default is "feature_summary".
#' @return A summary list with complex feature metrics.
#' @examples
#' ## Load example complex feature finding results:
#' complexFeatures <- exampleComplexFeatures
#' ## Run summary function:
#' summarizeFeatures(complexFeatures)
#' ## Filter complex features:
#' filteredComplexFeatures <- filterFeatures(complexFeatures)
#' ## Run summary function on filtered data:
#' summarizeFeatures(fileteredComplexFeatures)
#' @export
summarizeFeatures <- function(feature_table,plot=TRUE,PDF=FALSE,name="feature_summary"){
  features <- copy(feature_table)
  if ("protein_id" %in% names(features)){
    type <- "protein"
  } else if ("complex_id" %in% names(features)){
    type <- "complex"
  } else {
    stop("This is no protein or complex feature table. Please provide features from findComplexFeatures or findProteinFeatures.")
  }
  totalConfirmedHypotheses <- length(unique(features$complex_id))
  totalFeatures <- nrow(features)
  features[ , `:=`( COUNT = .N , IDX = 1:.N ) , by = complex_id ]
  setkey(features, "complex_id")
  feature_count_max <- unique(features)$COUNT
  summaryFeatureCount <- summary(feature_count_max) 
  totalHypothesesWithMultipleFeatures <- length(which(feature_count_max>=2))
  summaryCorrelation <- summary(features$peak_corr)
  summaryArea <- summary(features$area)
  if("mw_diff" %in% names(features)) {
    summaryMWdiff <- summary(features$mw_diff)
  }
  summaryNsubunitsAnnotated <- summary(features$n_subunits_annotated)
  summaryNsubunitsWithSignal <- summary(features$n_subunits_with_signal)
  summaryNsubunitsDetected <- summary(features$n_subunits_detected)
  summaryCompleteness <- summary(features$completeness)
  if("mw_diff" %in% names(features)) {
    summary <- list(type = type,
                   totalFeatures = totalFeatures,
                   totalConfirmedHypotheses = totalConfirmedHypotheses,
                   totalHypothesesWithMultipleFeatures = totalHypothesesWithMultipleFeatures,
                   summaryFeatureCount = summaryFeatureCount,
                   summaryCorrelation = summaryCorrelation,
                   summaryArea = summaryArea,
                   summaryMWdiff = summaryMWdiff,
                   summaryNsubunitsAnnotated = summaryNsubunitsAnnotated,
                   summaryNsubunitsWithSignal = summaryNsubunitsWithSignal,
                   summaryNsubunitsDetected = summaryNsubunitsDetected,
                   summaryCompleteness = summaryCompleteness
                   )
  } else {
    summary <- list(type = type,
                    totalFeatures = totalFeatures,
                    totalConfirmedHypotheses = totalConfirmedHypotheses,
                    totalHypothesesWithMultipleFeatures = totalHypothesesWithMultipleFeatures,
                    summaryFeatureCount = summaryFeatureCount,
                    summaryCorrelation = summaryCorrelation,
                    summaryArea = summaryArea,
                    summaryNsubunitsAnnotated = summaryNsubunitsAnnotated,
                    summaryNsubunitsWithSignal = summaryNsubunitsWithSignal,
                    summaryNsubunitsDetected = summaryNsubunitsDetected,
                    summaryCompleteness = summaryCompleteness
    )
  }
  if (plot) {
    features[,feature_id := .I]
    if("mw_diff" %in% names(features)) {
      plottingFeatures <- subset(features,select=c("feature_id","peak_corr","mw_diff","area","apex","n_subunits_detected","completeness")) 
    } else {
      plottingFeatures <- subset(features,select=c("feature_id","peak_corr","area","apex","n_subunits_detected","completeness")) 
    }
    plottingFeatures$area <- log(plottingFeatures$area)
    setnames(plottingFeatures,"area","log(area)")
    plottingFeatures[, names(plottingFeatures) := lapply(.SD, as.numeric)]
    plottingFeaturesMelt <- data.table::melt(plottingFeatures,id.vars=c("feature_id"))
    if (PDF) {
      pdf(gsub("$|\\.pdf$", ".pdf", name))
    }
    p <- ggplot(plottingFeaturesMelt, aes(x=value,fill=variable)) + 
      geom_histogram(bins = 50) + 
      facet_wrap(~variable, scales="free", ncol = 2) + 
      theme_bw() +
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + 
      guides(fill=FALSE) +
      ggtitle(paste0(type," feature summary")) + 
      theme(plot.title = element_text(size=14, face="bold"))
    print(p)
    plottingNsubcomplexes <- data.table(n_subcomplexes=feature_count_max)
    q <- ggplot(plottingNsubcomplexes,aes(x=n_subcomplexes)) +
      stat_bin(binwidth=1) +
      stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-0.5) +
      labs(x="N co-elution features",y="N complex hypotheses") +
      theme_classic() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      ggtitle(paste0(type," sub-feature summary")) + 
      theme(plot.title = element_text(size=14, face="bold"))
    print(q)
    if (PDF) {
      dev.off()
    }
  }
  return(summary)
}

