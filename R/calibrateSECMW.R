#' Calibrate SEC - establish connectivity between SEC fraction number and apparent
#' Molecular Weight (in kDa)
#'
#' @param std_weights_kDa The Molecular Weights of the standard proteins
#' @param std_elu_fractions The fraction numbers where these standard proteins eluted
#' @param Traces An Object of type Traces on which the calibration is to be performed
#' @return An additional column "apparentMW_kDa" in the Traces$fraction_annotation table
#' @export

#std_weights_kDa = c(1398, 699, 300, 150, 44, 17)
#std_elu_fractions = c(21, 31, 39, 48, 56.5, 63)

calibrateSECMW <- function(std_weights_kDa,std_elu_fractions,plot=TRUE) {
  calibrants <- NULL
  calibrants$logMW = log(std_weights_kDa)
  calibrants$fraction = std_elu_fractions
  if(length(calibrants$logMW) != length(calibrants$fraction)){
    stop("Check the input data. The vector of molecular weights and elution fractions need to be of same length.")
  }
  if(length(calibrants$logMW)<2 | length(calibrants$fraction)<2){
    stop("Check the input data. The vector of molecular weights and elution fractions need to be at least of length 2.")
  }
  # plot(std_weights_kDa, std_elu_fractions)
  model = lm(calibrants$logMW ~ calibrants$fraction, data = calibrants)
  # Visual check -> plot
  if(plot==TRUE){
    pdf("calibration.pdf")
      plot(calibrants$fraction, calibrants$logMW, xlab="SEC fraction", ylab="log10(MW)")
      lines(calibrants$fraction, fitted(model))
    dev.off()
  }
  # get linear model coefficients m and c for y = mx + c
  intercept <- as.numeric(model$coefficients[1])
  slope <- as.numeric(model$coefficients[2])

  # define resulting functions (how to export/update them for use by user?)
  MWtoSECfraction <- function(MW){
    round((log(MW)-intercept)/(slope), digits = 0)
  }
  SECfractionToMW <- function(SECfraction){
    exp(slope*SECfraction + intercept)
  }
  return(list(MWtoSECfraction=MWtoSECfraction,SECfractionToMW=SECfractionToMW))
}
