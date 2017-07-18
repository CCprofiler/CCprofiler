#' Complex feature detection
#' @description Run the sliding window algorithm to find complex features.
#' @param traces An object of class traces.
#' @param complex_hypothesis data.table containing complex hypotheses.
#' Should have the following columns:
#' \itemize{
#' \item complex_id: character strings, a unique id for every complex
#' \item complex_name: character strings, a unique name for every complex
#' \item protein_id: character strings, the protein id, e.g. Uniprot id
#' }
#' @param corr_cutoff Numeric, the correlation value for chromatograms above which
#'        peptides are considered to be coeluting, default=0.95.
#' @param window_size Integer, size of the window in fractions, default=15
#' @param parallelized Logical, if the computation should be done in parallel, default=\code{FALSE}.
#' @param n_cores Integer, the number of cores to use for parallel processing
#' (only applies if parallelized is \code{TRUE}), default=1.
#' @param collapse_method Character, method for collapsing multiple features into one feature:
#' \itemize{
#' \item "apex_only": collapses by apex
#' \item "apex_network": collapses by apex and connected network cluster}
#' Default="apex_only"
#' @param perturb_cutoff Numeric, the quantile to use in estimating the perturbation level, default="5%".
#'        Intensity values that are zero are replaced with random values that are
#'        below the specified quantile of the input values. Alternatively a
#'        cutoff value can be specified as an upper limit for perturbation values.
#'        This is nescessary for correlation calculation.
#' @param rt_height Numeric, RT cutoff for collapsing features. Defaults to 5.
#' @param smoothing_length Numeric, smoothing length of Savitzky-Golay filter. Defaults to 11.
#' @return A data.table containing protein complex features.
#' @export
#' @examples
#' ## Load example data
#' proteinTraces <- exampleProteinTraces
#' complexHypotheses <- exampleComplexHypotheses
#'
#' ## Perform co-elution signal detection
#' complexFeatures <- findComplexFeatures(traces=proteinTraces,
#'                                        complex_hypothesis=complexHypotheses)
#'
#' ## Inspect complex features
#' summarizeComplexFeatures(complexFeatures)
#'

findComplexFeatures <- function(traces,
                                complex_hypothesis,
                                corr_cutoff=0.95,
                                window_size=15,
                                parallelized=FALSE,
                                n_cores=1,
                                collapse_method="apex_only",
                                perturb_cutoff = "5%",
                                rt_height=5,
                                smoothing_length=11) { #MOD noise quantile can be user defined

  ## All complexes used for input
  inputComplexes <- unique(complex_hypothesis$complex_id)

  ## A helper function to execute the sliding window algorithm for a
  ## specific query complex.
  runSlidingWindow <- function(complex.id, traces, traces.imputed) {
    #complex.id <- inputComplexes[i]
    #cat(sprintf('CHECKING RUN:  %d / %d', i, length(inputComplexes)), '\n')
    # Extract the protein traces belonging to the current complex
    complex.subunits <- complex_hypothesis[complex_id == complex.id,
                                              protein_id]
    traces.subs <- subset(traces=traces,trace_subset_ids=complex.subunits)
    traces.imputed.subs <- traces.imputed[rownames(traces.imputed) %in% complex.subunits,,drop=F]
    complex.annotation <- subset(complex_hypothesis,complex_id == complex.id)
    # Run the algorithm
    try({
      if (nrow(traces.imputed.subs) >= 2) {
        complexFeaturesSW <- findComplexFeaturesSW(trace.mat=traces.imputed.subs,
                                                            corr.cutoff=corr_cutoff,
                                                            window.size=window_size)
        if((dim(complexFeaturesSW$features)[1] == 0) & (dim(complexFeaturesSW$features)[2] == 0)){
          return(list())
        }
        complexFeaturesPP <- findComplexFeaturesPP(traces.obj=traces.subs,
                                                   complexFeaturesSW=complexFeaturesSW,
                                                   smoothing_length=smoothing_length,
                                                   rt_height=rt_height)
        complexFeaturesCollapsed <- collapseComplexFeatures(complexFeature=complexFeaturesPP,rt_height=rt_height,collapse_method=collapse_method)
        if(dim(complexFeaturesCollapsed$features)[1] == 0){
          return(list())
        }
        # Calculate within peak boundary correlation
        complexFeaturesCollapsed.corr <- calculateFeatureCorrelation(traces.imputed.subs, complexFeaturesCollapsed)

        complexFeatureStoichiometries <- estimateComplexFeatureStoichiometry(traces.obj=traces.subs,
                                                                             complexFeaturesPP=complexFeaturesCollapsed.corr)
        complexFeatureAnnotated <- annotateComplexFeatures(traces,complexFeatureStoichiometries,complex.annotation)
        complexFeatureAnnotated
      } else {
        list()
      }
    })
  }


  #MOD calculate traces_matrix with imputed noise here and save it

  ## Impute noise for missing intensity measurements globally for all traces alternative
  trace.mat.imputed <- getIntensityMatrix(traces)
  n.zero.entries <- sum(trace.mat.imputed == 0) # number of ZERO values in matrix
  measure.vals <- trace.mat.imputed[trace.mat.imputed != 0]
  if(class(perturb_cutoff) == "character"){
    qt <- as.numeric(gsub("%","",perturb_cutoff))/100
    perturb_cutoff <- quantile(measure.vals, qt)
  }
  set.seed(123) # set seed to always get same results
  trace.mat.imputed[trace.mat.imputed == 0] <- sample(1:perturb_cutoff,size = n.zero.entries,
                                                      replace = TRUE)

  ##################################


  ## Execute the sliding window algorithm for each query complex.
  ## This computation can optionally be parstr(swf_ allelized.
  if (parallelized) {
    cl <- snow::makeCluster(n_cores)
    # setting a seed is absolutely crutial to ensure reproducible results!!!!!!!!!!!!!!!!!!!
    clusterSetRNGStream(cl,123)
    doSNOW::registerDoSNOW(cl)
    pb <- txtProgressBar(max = length(inputComplexes), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    sw.results <- foreach(i=seq_along(inputComplexes),
                          .packages=c('data.table', 'SECprofiler'),.options.snow = opts) %dopar% {
                            query.complex.id <- inputComplexes[i]
                            runSlidingWindow(query.complex.id, traces=traces,traces.imputed = trace.mat.imputed)
                          }
    close(pb)
    parallel::stopCluster(cl)
  } else {
    pb <- txtProgressBar(max = length(inputComplexes), style = 3)
    sw.results <- foreach(i=seq_along(inputComplexes)) %do% {
      setTxtProgressBar(pb, i)
      query.complex.id <- inputComplexes[i]
      #cat(sprintf('CHECKING RUN:  %d / %d', i, length(inputComplexes)), '\n')
      runSlidingWindow(query.complex.id, traces=traces, traces.imputed = trace.mat.imputed)
    }
    close(pb)
  }
  #names(sw.results) <- inputComplexes

  # @NEW remove results for complexes with ERROR message
  sel_errors <- which(sapply(sw.results, typeof) == "character")
  for (error_idx in sel_errors){
    sw.results[[error_idx]] <- list()
  }

  res <- list(sw.results=sw.results,
              inputComplexes=inputComplexes,
              corr.cutoff=corr_cutoff,
              window.size=window_size)
  class(res) <- 'complexFeatures'
  #res

  # only report table to make things easier
  complexRes <- .resultsToTable(res)
  complexRes

}

#' Reformat ComplexFeature object to data.table
#' @description Reformat complex feature object to a data.table.
#' @param complexFeatures A list containing various results.
#'         \itemize{
#'          \item \code{sw.results} A list of results of the function
#'                \code{findComplexFeatures}. One for each query complex.
#'          \item \code{input.complexes} A character vector of all query
#'                 complexes.
#'          \item \code{corr.cutoff} The correlation cutoff used.
#'          \item \code{window.size} The window size used.
#'          }
#' @return A data.table with all detected complex features of all query complexes
.resultsToTable <- function(complexFeatures){
  res <- complexFeatures$sw.results
  res.list <- lapply(seq(1:length(res)), function(i){
    features <- res[[i]]$features
    features
  })
  res <- do.call("rbind", res.list)
  res
}
