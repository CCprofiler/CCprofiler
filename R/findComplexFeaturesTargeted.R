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
                                        n.cores=parallel::detectCores()) {
    # Merge the association information into the trace DT to produce
    # a large DT where each row corresponds to the trace of a complex subunit.
    setkey(protein.traces, protein_id)
    setkey(complex.protein.assoc, protein_id)
    protein.traces.c <- merge(protein.traces,
                              complex.protein.assoc[, list(protein_id, complex_id)],
                              by='protein_id',
                              allow.cartesian=TRUE)
    # Generate decoys
    set.seed(42)
    protein.traces.decoy <- data.table::copy(protein.traces.c)
    protein.traces.decoy[, complex_id := sample(complex_id)]
    protein.traces.decoy[, complex_id := paste0('DECOY_', complex_id)]

    protein.traces.decoy[, is_decoy := TRUE]
    protein.traces.c[, is_decoy := FALSE]

    ## Concatenate target and decoy complexes
    protein.traces.all <- rbind(protein.traces.decoy, protein.traces.c)
    # Throw out all complexes with less than 2 traces
    protein.traces.all[, n_subunits_observed := length(protein_id), by='complex_id']
    protein.traces.all <- protein.traces.all[n_subunits_observed >= 2]
    protein.traces.all[, n_subunits_observed := NULL]

    ## All complexes used for input
    input.complexes.sw <- unique(protein.traces.all$complex_id)

    cl <- snow::makeCluster(n.cores)
    doSNOW::registerDoSNOW(cl)

    setkey(protein.traces.all, complex_id)

    sw.results <- foreach(i=seq_along(input.complexes.sw)) %do% {
        complex.id <- input.complexes.sw[i]
        cat('********************************************************************************\n')
        cat(sprintf('CHECKING RUN:  %d / %d', i, length(input.complexes.sw)), '\n')
        cat('********************************************************************************\n')

        # Extract the protein traces belonging to the current complex
        traces.subs <- protein.traces.all[complex_id == complex.id]
        is.decoy.complex <- traces.subs$is_decoy[1]
        protein.names <- traces.subs$protein_id
        # Convert traces to a matrix
        traces.mat <- as.matrix(subset(traces.subs,
                                       select=-c(complex_id, protein_id, is_decoy)))

        # Run the algorithm
        try({
            findComplexFeatures(traces.mat, protein.names,
                                corr.cutoff=corr.cutoff,
                                window.size=window.size)
        })
    }

    res <- list(sw_results=sw.results,
                input_complexes=input.complexes.sw,
                corr.cutoff=corr.cutoff,
                window.size=window.size,
                protein.complex.assoc=subset(protein.traces.all,
                                             select=c(protein_id, complex_id)))

    class(res) <- 'targetedComplexFeaturesSW'

    res
}
