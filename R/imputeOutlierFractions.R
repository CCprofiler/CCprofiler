#' Impute outlier fractions
#' @description Imputes values fo fractions that were labeled as outliers by getFractionInformation.
#' @param traces An object of type traces. These should be the raw, unprocessed traces.
#' @param impute_all logical, if all traces values should be imputed, default=FALSE
#' @param bound_left Numeric integer, the minimum number of non-zero values to the
#' left of a missing value to be replaced with \code{NA}.
#' @param bound_right Numeric integer, the minimum number of non-zero values to the
#' right of a missing value to be replaced with \code{NA}.
#' @param consider_borders Logical, whether to find missin values in the first/last fractions
#' (e.g. \code{0-1-1 -> NA-1-1} and \code{1-0-1-1 -> 1-NA-1-1} if the leftmost value is the
#' first fraction of the traces object).
#' @param method Character string, which interpolation method to use.
#' @details Unlike with standard imputation methods, missing values on the borders
#' are linearly extrapolated from the neighboring 2 values. Any imputed value below 0
#' is set to 0.
#' @return Numeric matrix with intensity values.
#' @examples
#' intensityMatrix <- getIntensityMatrix(examplePeptideTraces)
#' head(intensityMatrix)
#' @export
imputeOutlierFractions <- function(traces,
                              impute_all = FALSE,
                              bound_left = 2,
                              bound_right = 2,
                              consider_borders = TRUE,
                              method = c("mean", "spline")){
  UseMethod("imputeOutlierFractions", traces)
}

#' @describeIn imputeOutlierFractions Impute outlier fractions in traces object
#' @export
imputeOutlierFractions.traces <- function(traces,
                                          impute_all = FALSE,
                                          bound_left = 2,
                                          bound_right = 2,
                                          consider_borders = TRUE,
                                          method = c("mean", "spline"), ...){
  .tracesTest(traces)

  # get outlier ids
  outlier_fractions <- traces$fraction_annotation[isOutlier==TRUE]$id
  # remove values from outlier fractions
  flagedTraces <- copy(traces)
  if (impute_all==TRUE) {
    flagedTraces$traces[,(outlier_fractions):= as.numeric(NA),by=id]
  } else {
    flagedTraces$traces[,(outlier_fractions):=0,by=id]
    flagedTracesToImpute <- findMissingValues(flagedTraces)
  }
  # run imputation function
  resTraces <- imputeMissingVals(flagedTracesToImpute,method="spline")

  .tracesTest(resTraces)
  return(resTraces)
}

#' @describeIn imputeOutlierFractions Impute outlier fractions in multiple traces objects
#' @export
imputeOutlierFractions.tracesList <- function(tracesList,
                                              impute_all = FALSE,
                                              bound_left = 2,
                                              bound_right = 2,
                                              consider_borders = TRUE,
                                              method = c("mean", "spline"), ...){
  .tracesListTest(tracesList)
  tracesListRes <- lapply(tracesList, imputeOutlierFractions.traces,
                          impute_all = FALSE,
                          bound_left = 2,
                          bound_right = 2,
                          consider_borders = TRUE,
                          method = c("mean", "spline"), ...)
  class(tracesListRes) <- "tracesList"
  .tracesListTest(tracesListRes)
  return(tracesListRes)
}
