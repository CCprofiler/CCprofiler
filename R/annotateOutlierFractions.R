
#' Annotate outlier fractions
#' @param traces Object of class traces.
#' @param excludeBoarder If boarder fractions should be excluded as outliers, default is TRUE.
#' @return Object of class traces with annotation of fractions and outliers.
#' @export
#' @example
#' inputTraces <- examplePeptideTraces
#' exampleTracesAnnotatedOutliers <- annotateOutlierFractions(inputTraces)
#' exampleTracesAnnotatedOutliers$fraction_annotation
annotateOutlierFractions <- function(traces, excludeBoarder=TRUE,
                          plot = TRUE, PDF=FALSE, name="OutlierFractions"){
  UseMethod("annotateOutlierFractions", traces)
}

#' @describeIn annotateOutlierFractions Detect outlier fractions in single traces object
#' @export
annotateOutlierFractions.traces <- function(traces, excludeBoarder=TRUE,
                          plot = TRUE, PDF=FALSE, name="OutlierFractions", ...){
  .tracesTest(traces)
  traces_res <- copy(traces)
  traces_updated <- updateTraces(traces_res)
  dist_data <- subset(traces_updated$fraction_annotation, select=c("missingValues","loessResidualsIdCounts","loessResiduals"))
  dfN <- 3
  if (excludeBoarder) {
    dist_data$missingValues[1] = 0
    dist_data$missingValues[length(dist_data$missingValues)] = 0
  }
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  if (max(dist_data$missingValues)==0) {
    dist_data[,missingValues := NULL]
    dfN <- 2
  } else {
    dist_data[, missingValues := range01(missingValues)]
  }
  dist_data[, loessResiduals := abs(loessResiduals)]
  dist_data[, loessResiduals := range01(loessResiduals)]
  dist_data[, loessResidualsIdCounts := abs(loessResidualsIdCounts)]
  dist_data[, loessResidualsIdCounts := range01(loessResidualsIdCounts)]
  mahalanobis_dist <- mahalanobis(dist_data, colMeans(dist_data), cov(dist_data), tol = 1e-17)
  traces_updated$fraction_annotation[,mahalanobis := mahalanobis_dist]
  traces_updated$fraction_annotation[,chisq.pvalue := pchisq(mahalanobis_dist, df=dfN, lower.tail=FALSE)]
  traces_updated$fraction_annotation[,isOutlier := ifelse(chisq.pvalue < 0.01, TRUE, FALSE)]
  if (excludeBoarder==TRUE) {
    traces_updated$fraction_annotation$isOutlier[1] = FALSE
    traces_updated$fraction_annotation$isOutlier[length(traces_updated$fraction_annotation$isOutlier)] = FALSE
  }
  if (plot) {
    message("plot")
    if (PDF) {
      pdf(paste0(name,".pdf"))
    }
    chi3sq <- sort(qchisq(ppoints(nrow(traces_updated$fraction_annotation)), df = dfN))
    D2 <- sort(traces_updated$fraction_annotation$mahalanobis)
    dt <- data.table(chi3sq=chi3sq, D2=D2)
    qqplot <- ggplot(dt, aes(x=chi3sq,y=D2)) + geom_point() +
        geom_abline() + ggtitle("Q-Q plot of Mahalanobis")
    print(qqplot)
    d2plot <- ggplot(traces_updated$fraction_annotation, aes(x=id, y=mahalanobis)) +
      geom_point() + ggtitle("Mahalanobis")
    print(d2plot)
    d_data <- dist_data
    d_data[,id := .I]
    if ("missingValues" %in% names(dist_data)){
      d1 <- ggplot(dist_data, aes(x=id, y=missingValues)) +
        geom_point() + ggtitle("missingValues")
      print(d1)
    }
    d2 <- ggplot(dist_data, aes(x=id, y=loessResidualsIdCounts)) +
      geom_point() + ggtitle("loessResidualsIdCounts")
    print(d2)
    d3 <- ggplot(dist_data, aes(x=id, y=loessResiduals)) +
      geom_point() + ggtitle("loessResiduals")
    print(d3)
    if(PDF){
      dev.off()
    }
  }
  .tracesTest(traces_updated)
  return(traces_updated)
}

#' @describeIn annotateOutlierFractions Detect outlier fractions in multiple traces objects
#' @export

annotateOutlierFractions.tracesList <- function(traces, excludeBoarder=TRUE,
                          plot = TRUE, PDF=FALSE, name="OutlierFractions", ...){
  .tracesListTest(traces)
  res <- lapply(names(traces), function(sample){
          message(paste("Annotate outliers for sample", sample, "..."))
          traces_i <- traces[[sample]]
          name <- paste0(name, "_", sample)
          annotateOutlierFractions.traces(traces = traces_i, excludeBoarder=excludeBoarder,
                          plot = plot, PDF=PDF, name=name, ...)
          })
  names(res) <- names(traces)
  class(res) <- "tracesList"
  .tracesListTest(res)
  return(res)
}
