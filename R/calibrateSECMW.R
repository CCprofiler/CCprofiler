#' calibrateSECMW
#' @description Calibrate SEC - establish connectivity between SEC fraction number and apparent
#' Molecular Weight (in kDa)
#' @param calibration_table tab-separated calibration table, file or R data.table
#' Columns: 
#' std_weights_kDa = The Molecular Weights of the standard proteins
#' std_elu_fractions = The fraction numbers where these standard proteins eluted
#' @param plot logical, if to plot calibration
#' @param PDF logical, if to produce a PDF
#' @return List of two functions: MWtoSECfraction and SECfractionToMW
#' @export

#std_weights_kDa = c(1398, 699, 300, 150, 44, 17)
#std_elu_fractions = c(19, 29, 37, 46, 54.5, 61)

calibrateSECMW <- function(calibration_table,plot=TRUE,PDF=FALSE) {
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
  model = lm(calibrants$logMW ~ calibrants$fraction, data = calibrants)
  # Visual check -> plot
  if(plot){
    if(PDF){
      pdf("calibration.pdf")
    }
    plot(calibrants$fraction, calibrants$logMW, xlab="SEC fraction", ylab="log10(MW)",main="calibration")
    lines(calibrants$fraction, fitted(model))
    if(PDF){
      dev.off()
    }
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
