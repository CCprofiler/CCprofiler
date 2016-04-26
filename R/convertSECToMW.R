#' Convert a fraction number to an estimated molecular weight.
#'
#' @param sec The fraction number.
#' @return A molecular weight estimate.
#' @export
convertSECToMW <- function(sec) {
    logMW <- -0.0453 * sec + 4.205
    exp(logMW)
}
