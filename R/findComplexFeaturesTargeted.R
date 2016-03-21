#' Run the sliding window algorithm for a number of protein complexes.
#' 
#' @param protein.traces A wide format data.table where the first column
#'        is named \code{protein_id} and all other columns correspond to
#'        intensity values measured in SEC fractions. 
#' @param complex.protein.assoc A data.table that encodes the input complexes
#'        and their subunit compositions. The data.table should have the
#'        format:
#'        \itemize{
#'         \item \code{complex_id} A character vector that uniquely identifies
#'               a complex.
#'         \item \code{protein_id} A character vector that uniquely identifies
#'               a protein and links it to a trace contained in
#'               \code{protein.traces}.
#'        }
#' @param corr.cutoff The correlation cutoff to use in the sliding window
#'        algorithm. Optional.
#' @param window.size The window size to use in the sliding window algorithm.
#'        Optional.
#' @param parallelized If the computation should be done in parallel.
#'        Optional.
#' @param n.cores The number of cores to use for parallel processing. Optional.
#' @return A list containing various results.
#'         \itemize{
#'          \item \code{sw.results} A list of results of the function
#'                \code{findComplexFeatures}. One for each query complex.
#'          \item \code{complex.stats} A data.table specifying for each query
#'                complex how many of its subunits were detected.
#'          \item \code{input.complexes} All query complexes.
#'          \item \code{input.complexes} A character vector of all query
#'                 complexes.
#'          \item \code{corr.cutoff} The correlation cutoff used.
#'          \item \code{window.size} The window size used.
#'          \item \code{complex.protein.assoc} The protein-complex associations
#'                used as the original input argument.
#'         }
#' @examples
#' # NOT RUN:
#' # sample.complexes <-
#' #     sample(unique(corum.complex.protein.assoc$complex_id), 10)
#' # complex.protein.assoc <-
#' #     corum.complex.protein.assoc[complex_id %in% sample.complexes]
#' # findComplexFeaturesTargeted(protein.traces, complex.protein.assoc)
#' @export
findComplexFeaturesTargeted <- function(protein.traces,
                                        complex.protein.assoc,
                                        corr.cutoff=0.99,
                                        window.size=15,
                                        parallelized=FALSE,
                                        n.cores=parallel::detectCores()) {
    setkey(protein.traces, protein_id)
    setkey(complex.protein.assoc, complex_id)

    ## All complexes used for input
    input.complexes <- unique(complex.protein.assoc$complex_id)

    ## A helper function to execute the sliding window algorithm for a 
    ## specific query complex.
    runSlidingWindow <- function(complex.id) {
        complex.id <- input.complexes[i]
        cat(sprintf('CHECKING RUN:  %d / %d', i, length(input.complexes)), '\n')
        # Extract the protein traces belonging to the current complex
        complex.subunits <- complex.protein.assoc[complex_id == complex.id,
                                                  protein_id]
        traces.subs <- protein.traces[protein_id %in% complex.subunits]
        protein.names <- traces.subs$protein_id
        # Convert traces to a matrix
        traces.mat <- as.matrix(subset(traces.subs, select=-protein_id))
        # Run the algorithm
        try({
            if (nrow(traces.mat) >= 2) {
                findComplexFeatures(traces.mat, protein.names,
                                    corr.cutoff=corr.cutoff,
                                    window.size=window.size)
            } else {
                list()
            }
        })
    }

    ## Execute the sliding window algorithm for each query complex.
    ## This computation can optionally be parallelized.
    if (parallelized) {
        cl <- snow::makeCluster(n.cores)
        doSNOW::registerDoSNOW(cl)
        sw.results <- foreach(i=seq_along(input.complexes),
                             .packages=c('data.table', 'SECprofiler')) %dopar% {
            query.complex.id <- input.complexes[i]
            runSlidingWindow(query.complex.id)
        }
    } else {
        sw.results <- foreach(i=seq_along(input.complexes)) %do% {
            query.complex.id <- input.complexes[i]
            runSlidingWindow(query.complex.id)
        }
    }
    names(sw.results) <- input.complexes

    ## Calculate complex detection statistics such as how many subunits
    ## were detected in a subgroup.
    complex.stats <- complex.protein.assoc[, list(n_subunits_annotated=.N),
                                           by=complex_id]
    complex.stats[, n_subunits_detected := 0]
    setkey(complex.stats, complex_id)
    for (i in seq_along(input.complexes)) {
        complex.id <- input.complexes[i]
        res <- sw.results[[i]]
        if (!is.null(res$subgroups.dt) && nrow(res$subgroups.dt) > 0) {
            max.detected.subunits <- max(res$subgroups.dt$n_subunits)
            complex.stats[complex_id == as.character(complex.id),
                          n_subunits_detected := max.detected.subunits]
        }
    }
    complex.stats[, completeness := n_subunits_detected / n_subunits_annotated]

    res <- list(sw.results=sw.results,
                complex.stats=complex.stats,
                input.complexes=input.complexes,
                corr.cutoff=corr.cutoff,
                window.size=window.size,
                complex.protein.assoc=complex.protein.assoc)

    class(res) <- 'targetedComplexFeaturesSW'

    res
}
