
#' calibrate Molecular Weight
#' @description Establish connectivity between fraction number and apparent
#' Molecular Weight(MW) (in kDa). Many fractionation based PCP expeiments will have a set
#' of standard proteins with known MW spiked in the samples. These reference
#' proteins can be used in this function to calibrate the corresponding MW of every fraction.
#' While this is not neccesary for downstream analysis, it is recommended if standard proteins
#' are available, since it facilitates the process of drawing biological conclusions.
#' @details The function uses a standard \code{lm()} of the form \code{log2(MW) ~ fraction}.
#' Therefore, the fractionation technique is required to have a log-linear relationship between
#' fraction and molecular weight (True for e.g. standard SEC methods).
#' The estimated slope and intersect of the model are then used for the returned calibration functions.
#' @param calibration_table A table with protein standards, file or R data.table.
#' Columns:
#' \itemize{
#'  \item \code{std_weights_kDa} Numeric, the MW of the standard proteins.
#'  \item \code{std_elu_fractions} Numeric, the fraction numbers where these standard proteins eluted.
#' }
#' For an example see \code{exampleCalibrationTable}.
#' @param plot Logical, whether to plot calibration curve. Defaults to \code{TRUE}.
#' @param PDF Logical, whether to produce a PDF file in the working directory.
#' Defaults to \code{FALSE}.
#' @return List of two functions: MWtoFraction and FractionToMW.
#' @export
#' @examples 
#'
#' ## Load example data
#' standardProteinTable <- exampleCalibrationTable
#' ## Regress the standard proteins against their elution fractions
#' calibrationFunctions <- calibrateMW(standardProteinTable,
#'                                        plot = TRUE,
#'                                        PDF = FALSE)
#' # The provided plot can be used for quality assessment. The points represent the
#' # standard proteins and should agree with the model (solid line))
#'
#' ## The molecular weight (in kDa) of an arbitrary fraction can now be calculated.
#' calibrationFunctions$FractionToMW(14)
#'
#' ## Arbitrary MWs (in kDa) can also be converted to fraction numbers
#' calibrationFunctions$MWtoFraction(3020) # i.e. 3020 kDa
#'

calibrateMW <- function(calibration_table,
                           plot=TRUE,
                           PDF=FALSE) {
  # use or if path to file read annotation table and add id column
  if (class(calibration_table)[1] == "character") {
    if (file.exists(calibration_table)) {
      message('reading calibration table ...')
      calibration_table  <- data.table::fread(calibration_table, header = TRUE)
    } else {
      stop("calibration_table file doesn't exist")
    }
  } else if (all(class(calibration_table) != c("data.table","data.frame"))) {
    stop("calibration_table input is neither file name or data.table")
  }
  if(nrow(calibration_table) == 0) stop("calibration_table seems to be empty")
  calibrants <- NULL
  calibrants$logMW = log(calibration_table$std_weights_kDa)
  calibrants$fraction = calibration_table$std_elu_fractions
  if(length(calibrants$logMW) != length(calibrants$fraction)){
    stop("Check the input data. The vector of molecular weights and elution fractions need to be of same length.")
  }
  if(length(calibrants$logMW)<2 | length(calibrants$fraction)<2){
    stop("Check the input data. The vector of molecular weights and elution fractions need to be at least of length 2.")
  }
  # plot(std_weights_kDa, std_elu_fractions)
  modelMW = lm(calibrants$logMW ~ calibrants$fraction, data = calibrants)
  # Visual check -> plot
  if(plot){
    if(PDF){
      pdf("calibration.pdf")
    }
    plot(calibrants$fraction, calibrants$logMW, xlab="Fraction",
         ylab="log10(MW (kDa))", main="calibration")
    legend(x = min(calibrants$fraction),
           y = min(calibrants$logMW),
           legend = paste("R^2 = ", round(summary(modelMW)$r.squared,2)),
           yjust = 0)
    lines(calibrants$fraction, fitted(modelMW))
    if(PDF){
      dev.off()
    }
  }
  # get linear model coefficients m and c for y = mx + c
  intercept <- as.numeric(modelMW$coefficients[1])
  slope <- as.numeric(modelMW$coefficients[2])

  # define resulting functions (how to export/update them for use by user?)
  MWtoFraction <- function(MW){
    round((log(MW)-intercept)/(slope), digits = 0)
  }
  FractionToMW <- function(fraction){
    exp(slope*fraction + intercept)
  }
  return(list(MWtoFraction=MWtoFraction,FractionToMW=FractionToMW))
}
