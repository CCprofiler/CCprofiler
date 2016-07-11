#' Convert a fraction number to an estimated molecular weight.
#'
#' @param SECfraction The fraction number.
#' @return Apparent MW at this SEC fraction in kDa.
#' @export
convertSECToMW <- function(SECfraction) {
    logMW <- -0.1043329 * SECfraction + 9.682387
    exp(logMW)
}

