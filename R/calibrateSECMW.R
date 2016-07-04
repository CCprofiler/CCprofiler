#' Calibrate SEC - establish connectivity between SEC fraction number and apparent
#' Molecular Weight (in kDa)
#'
#' @param std_weights_kDa The Molecular Weights of the standard proteins
#' @param std_elu_fractions The fraction numbers where these standard proteins eluted
#' @param Traces An Object of type Traces on which the calibration is to be performed
#' @return An additional column "apparentMW_kDa" in the Traces$fraction_annotation table
#' @export

calibrateSECMW <- function(Traces,
                           std_weights_kDa = c(1398, 699, 300, 150, 44, 17),
                           std_elu_fractions = c(21, 31, 39, 48, 56.5, 63)) {
  calibrants <- NULL
  calibrants$logMW = log(std_weights_kDa)
  calibrants$fraction = std_elu_fractions
  # plot(std_weights_kDa, std_elu_fractions)
  model = lm(calibrants$logMW ~ calibrants$fraction, data = calibrants)
  # Visual check -> plot
  plot(calibrants$fraction, calibrants$logMW)
  lines(calibrants$fraction, fitted(model))
  # get linear model coefficients m and c for y = mx + c
  intercept <- model$coefficients[1]
  slope <- model$coefficients[2]
  
  # define resulting functions (how to export/update them for use by user?)
  MWtoSECfraction <- function(MW){
    round((log(MW)-intercept)/(slope), digits = 0)
  }
  SECfractionToMW <- function(SECfraction){
    exp(slope*SECfraction + intercept)
  }
  
  Traces$fraction_annotation[, apparentMW_kDa:= SECfractionToMW(id)]
  Traces
}

