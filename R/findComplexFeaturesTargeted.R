#' Run the sliding window algorithm over a wide format data.table of protein
#' traces.
#' @export
findComplexFeaturesTargeted <- function(protein.traces.sw,
                                        corr.cutoff=0.99,
                                        window.size=15,
                                        n.cores=detectCores()) {

    setkey(protein.traces.sw, complex_id)

    # Generate decoys
    set.seed(42)
    protein.traces.sw.decoy <- data.table::copy(protein.traces.sw)
    protein.traces.sw.decoy[, complex_id := sample(complex_id)]
    protein.traces.sw.decoy[, complex_id := paste0('DECOY_', complex_id)]
    protein.traces.sw.decoy[, is_decoy := TRUE]
    protein.traces.sw[, is_decoy := FALSE]

    ## Concatenate target and decoy complexes
    protein.traces.sw.all <- rbind(protein.traces.sw.decoy, protein.traces.sw)
    # Throw out all complexes with less than 2 traces
    protein.traces.sw.all[, n_subunits_observed := length(protein_id), by='complex_id']
    protein.traces.sw.all <- protein.traces.sw.all[n_subunits_observed >= 2]
    protein.traces.sw.all[, n_subunits_observed := NULL]

    ## All complexes used for input
    input.complexes.sw <- unique(protein.traces.sw.all$complex_id)

    registerDoMC(n.cores)

    sw.results <- foreach(i=seq_along(input.complexes.sw)) %dopar% {
        complex.id <- input.complexes.sw[i]
        cat('********************************************************************************\n')
        cat(sprintf('CHECKING RUN:  %d / %d', i, length(input.complexes.sw)), '\n')
        cat('********************************************************************************\n')

        # Extract the protein traces belonging to the current complex
        traces.subs <- protein.traces.sw.all[complex_id == complex.id]
        is.decoy.complex <- traces.subs$is_decoy[1]
        protein.names <- traces.subs$protein_id
        # Convert traces to a matrix
        traces.mat <- as.matrix(subset(traces.subs, select=-c(complex_id, protein_id, is_decoy)))

        # Run the algorithm
        try({
            detectSubgroupsSW(traces.mat, protein.names,
                              corr.cutoff=corr.cutoff,
                              window.size=window.size)
        })
    }

    res <- list(sw_results=sw.results,
                input_complexes=input.complexes.sw,
                corr.cutoff=corr.cutoff,
                window.size=window.size,
                protein.complex.assoc=subset(protein.traces.sw.all, select=c(protein_id, complex_id)
                ))

    class(res) <- 'targetedComplexFeaturesSW'

    res
}
