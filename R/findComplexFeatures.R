#' Run the sliding window algorithm for a number of protein complexes.
#'
#' @param traces.obj An object of type \code{traces.obj}.
#' @param complex.protein.assoc A data.table that encodes the input complexes
#'        and their subunit compositions. The data.table should have the
#'        format:
#'        \itemize{
#'         \item \code{complex_id} A character vector that uniquely identifies
#'               a complex.
#'         \item \code{complex_name} A character vector that contains
#'               the complex name.
#'         \item \code{protein_id} A character vector that uniquely identifies
#'               a protein and links it to a trace contained in
#'               \code{traces.obj}.
#'        }
#' @param corr.cutoff The correlation value for chromatograms above which
#'        proteins are considered to be coeluting.
#' @param window.size Size of the window. Numeric.
#' @param parallelized If the computation should be done in parallel.
#'        Optional.
#' @param n.cores The number of cores to use for parallel processing. Optional.
#' @return A list containing various results.
#'         \itemize{
#'          \item \code{sw.results} A list of results of the function
#'                \code{findComplexFeatures}. One for each query complex.
#'          \item \code{input.complexes} A character vector of all query
#'                 complexes.
#'          \item \code{corr.cutoff} The correlation cutoff used.
#'          \item \code{window.size} The window size used.
#'         }
#' @examples
#' # NOT RUN:
#' # protein.ids <- corum.complex.protein.assoc[complex_id == 181, protein_id]
#' # traces <- subset(protein.traces[protein_id %in% protein.ids],
#' #                  select=-protein_id)
#' # sample.complexes <- c(181, 191)
#' # complex.protein.assoc <-
#' #     corum.complex.protein.assoc[complex_id %in% sample.complexes]
#' # findComplexFeaturesSWBulk(traces,
#' #                           protein.ids,
#' #                           protein.mw.conc[protein_id %in% protein.ids],
#' #                           complex.protein.assoc)
#' @export
findComplexFeatures <- function(traces.obj,
                                  complex.protein.assoc,
                                  MWSECcalibrationFunctions,
                                  corr.cutoff=0.95,
                                  window.size=15,
                                  parallelized=FALSE,
                                  n.cores=parallel::detectCores()) {

    ## All complexes used for input
    input.complexes <- unique(complex.protein.assoc$complex_id)

    ## A helper function to execute the sliding window algorithm for a
    ## specific query complex.
    runSlidingWindow <- function(complex.id) {
        #complex.id <- input.complexes[i]
        #cat(sprintf('CHECKING RUN:  %d / %d', i, length(input.complexes)), '\n')
        # Extract the protein traces belonging to the current complex
        complex.subunits <- complex.protein.assoc[complex_id == complex.id,
                                                  protein_id]
        traces.subs <- subset(traces.obj=traces.obj,trace_ids=complex.subunits)
        complex.annotation <- subset(complex.protein.assoc,complex_id == complex.id)
        # Run the algorithm
        try({
            if (nrow(traces.subs$traces) >= 2) {
              complexFeaturesSW <- findComplexFeaturesSW(traces.obj=traces.subs,
                                                          corr.cutoff=corr.cutoff,
                                                          window.size=window.size)
              if((dim(complexFeaturesSW$features)[1] == 0) & (dim(complexFeaturesSW$features)[2] == 0)){
                return(list())
              }
              complexFeaturesPP <- findComplexFeaturesPP(traces.obj=traces.subs,
                                                          complexFeaturesSW=complexFeaturesSW)
              complexFeatureStoichiometries <- estimateComplexFeatureStoichiometry(traces.obj=traces.subs,
                                                          complexFeaturesPP=complexFeaturesPP)
              complexFeatureAnnotated <- annotateComplexFeatures(traces.obj,complexFeatureStoichiometries,complex.annotation,MWSECcalibrationFunctions)
              complexFeatureAnnotated
            } else {
                list()
            }
        })
    }

    ## Execute the sliding window algorithm for each query complex.
    ## This computation can optionally be parstr(swf_ allelized.
    if (parallelized) {
      cl <- snow::makeCluster(n.cores,outfile="")
      # setting a seed is absolutely crutial to ensure reproducible results!!!!!!!!!!!!!!!!!!!
      clusterSetRNGStream(cl,123)
      doSNOW::registerDoSNOW(cl)
      pb <- txtProgressBar(max = length(input.complexes), style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      sw.results <- foreach(i=seq_along(input.complexes),
        .packages=c('data.table', 'SECprofiler'),.options.snow = opts) %dopar% {
        query.complex.id <- input.complexes[i]
        runSlidingWindow(query.complex.id)
      }
      close(pb)
      parallel::stopCluster(cl)
    } else {
      pb <- txtProgressBar(max = length(input.complexes), style = 3)
      sw.results <- foreach(i=seq_along(input.complexes)) %do% {
        setTxtProgressBar(pb, i)
        query.complex.id <- input.complexes[i]
        #cat(sprintf('CHECKING RUN:  %d / %d', i, length(input.complexes)), '\n')
        runSlidingWindow(query.complex.id)
      }
      close(pb)
    }
    #names(sw.results) <- input.complexes

    # @NEW remove results for complexes with ERROR message
    sel_errors <- which(sapply(sw.results, typeof) == "character")
    for (error_idx in sel_errors){
      sw.results[[error_idx]] <- list()
    }

    res <- list(sw.results=sw.results,
                input.complexes=input.complexes,
                corr.cutoff=corr.cutoff,
                window.size=window.size)

    class(res) <- 'complexFeatures'

    res
}
