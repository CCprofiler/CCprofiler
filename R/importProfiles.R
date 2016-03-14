#' Import peptide profiles from an OpenSWATH experiment.
#'
#' @param file.name Quantitative MS data result file.
#' @param annotation.table A data frame with two columns `run_id` and
#'        `sec_fraction` that map the run identifier as contained in `file.name`
#'        to a SEC elution fraction.
#' @return An object of class SECprofiles that can be processed with the herein
#'         contained function.s
#'
#' @export
importProfilesFromOpenSWATH <- function(file.name, annotation.table) {
    # TODO: Don't load whole table into RAM since table contains other
    # information besides the intensity measurements.
    # ...
    traces <- list()
    class(traces) <- 'SECprofiles'
    traces
}

