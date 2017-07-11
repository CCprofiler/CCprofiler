#' Summarize co-elution features
#' @description Summarize co-elution features
#' @param feature_table data.table as reported by \code{\link[SECprofiler]{findComplexFeatures}} or \code{\link[SECprofiler]{findProteinFeatures}}.
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
summarizeFeatures <- function(feature_table,plot=TRUE,pdf=FALSE,name="feature_summary"){
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
  summaryMWdiff <- summary(features$mw_diff)
  summaryNsubunitsAnnotated <- summary(features$n_subunits_annotated)
  summaryNsubunitsWithSignal <- summary(features$n_subunits_with_signal)
  summaryNsubunitsDetected <- summary(features$n_subunits_detected)
  summaryCompleteness <- summary(features$completeness)
  list(type = type,
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
  if (plot) {
    features[,feature_id := .I]
    plottingFeatures <- subset(features,select=c("feature_id","peak_corr","area","mw_diff","n_subunits_detected","completeness")) 
    plottingFeatures$area <- log(plottingFeatures$area)
    plottingFeaturesMelt <- data.table::melt(plottingFeatures,id.vars=c("feature_id"))
    if (PDF) {
      pdf(gsub("$|\\.csv$", ".csv", name))
    }
    p <- ggplot(plottingFeaturesMelt, aes(x=value)) + geom_histogram(bins = 50) + facet_wrap(~variable, scales="free", ncol = 2) 
    print(p)
    if (PDF) {
      dev.off()
    }
  }
}
