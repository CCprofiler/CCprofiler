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
#' @return A list containing the individual results of the sliding window
#'         algorithm, one for each complex contained in
#'         \code{complex.protein.assoc}.
#' @examples
#' # NOT RUN:
#' # findComplexFeaturesTargeted(protein.traces, corum.complex.protein.assoc)
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

    runSlidingWindow <- function(i) {
        complex.id <- input.complexes[i]
        cat('********************************************************************************\n')
        cat(sprintf('CHECKING RUN:  %d / %d', i, length(input.complexes)), '\n')
        cat('********************************************************************************\n')
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

    if (parallelized) {
        cl <- snow::makeCluster(n.cores)
        doSNOW::registerDoSNOW(cl)
        sw.results <- foreach(i=seq_along(input.complexes),
                              .export='data.table') %dopar% {
            runSlidingWindow(i);
        }
    } else {
        sw.results <- foreach(i=seq_along(input.complexes)) %do% {
            runSlidingWindow(i);
        }
    }

    names(sw.results) <- input.complexes

    res <- list(sw.results=sw.results,
                input.complexes=input.complexes,
                corr.cutoff=corr.cutoff,
                window.size=window.size,
                complex.protein.assoc=complex.protein.assoc)

    class(res) <- 'targetedComplexFeaturesSW'

    res
}
