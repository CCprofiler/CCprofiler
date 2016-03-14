#' Detect subgroups within a window. This is a helper function and should 
#' be called by `findComplexFeatures`.
#'
#' @param sec.profiles An object of class SECprofiles that is to be filtered.
#' @param stretch.len The minimal length a stretch of contiguous identifications
#'        has to have in order not to be removed.
#' @return An object of class SECprofiles that represents the filtered
#'         chromatograms.
#'
#' @export
filterIdStretches <- function(sec.profiles, stretch.len) {
    # ...
}

