
#' A helper function to extend a list of complex features with additional
#' information.
#'
#' @param features A data.table of complex feature candidates with the
#'        following format:
#'        \itemize{
#'         \item \code{subgroup} A semicolon-separated list of protein
#'                               identifiers.
#'         \item \code{left_sw} The left boundary of the feature.
#'         \item \code{right_sw} The right boundary of the feature.
#'         \item \code{score} The intra-feature correlation.
#'        }
#' @param trace.mat A matrix where rows correspond to protein traces.
#'        This is the matrix that was used to find complex features.
#' @param protein.names A character vector specifying the protein identifiers
#'        belonging to the rows in \code{trace.mat}.
#' @param protein.mw.conc A data.table that stores the molecular weight and
#'        estimate of the absolute abundance for each subunit.
#'        \itemize{
#'         \item \code{protein_id}
#'         \item \code{protein_mw}
#'         \item \code{protein_concentration}
#'        }
#' @param protein.abundances A numeric vector holding estimates of the absolute
#'        abundance of each subunit.
#' @return The same data.table as the input argument extended with the
#'         following columns:
#'         \itemize{
#'          \item \code{n_subunits} The number of subunits in the feature.
#'          \item \code{stoichiometry} The intensity-based stoichiometry.
#'          \item \code{mw_estimated} The estimated molecular weight.
#'          \item \code{mw_apparent} The apparent mw.
#'          \item \code{mw_delta} The delta mw.
#'         }
extendComplexFeatures <- function(features, trace.mat,
                                  protein.names,
                                  protein.mw.conc) {
    # Compute the number of subunits in each complex feature
    features[, n_subunits := length(strsplit(as.character(subgroup), ';')[[1]]),
             by=subgroup]
    # Estimate the molecular weight of the complex by using information
    # gathered by calibrating the column with molecules of known molecular
    # weight. This is known as the 'apparent' molecular weight.
    features[, mw_apparent := convertSECToMW(left_sw+((right_sw - left_sw)/2)),
             by=subgroup] # @CORRECTED for central sec fraction

    # Produce a long list version of the trace matrix since its more convenient
    # for downstream processing.
    traces.dt <- data.table(trace.mat)
    traces.dt[, protein_id := protein.names]
    traces.long <- melt(traces.dt, id.vars='protein_id', value.name='intensity',
                        variable.name='fraction', variable.factor=FALSE)
    traces.long[, fraction := as.numeric(fraction)]
    setkey(traces.long, protein_id)

    # Next, create a data.table that holds information on subunits.
    protein.info <- traces.long[, list(total_intensity=sum(intensity)), by=protein_id]
    # Merge the data.table holding the molecular weight and abs. abundance
    # estimates and the data.table with the total intensity of each protein.
    setkey(protein.mw.conc, protein_id)
    setkey(protein.info, protein_id)
    protein.info <- protein.info[protein.mw.conc]
    setkey(protein.info, protein_id)

    # Loop over each feature and try to estimate the stoichiometry of its
    # complex as well as its molecular mass.
#    for (i in 1:nrow(features)) {
#        feature <- features[i]
#        subunits <- strsplit(feature$subgroup, ';')[[1]]
#        # Extract only the traces for the subunits that make up this feature.
#        subunit.traces <- traces.long[subunits]
#        # Extract the relevant information from the protein.info DT.
#        subunit.mws <- protein.info[subunits, protein_mw]
#        subunit.abundances <- protein.info[subunits, protein_concentration]
#        subunit.total.intensities <- protein.info[subunits, total_intensity]
#        # complex and its stoichiometry.
#        res <- estimateComplexMass(subunit.traces,
#                                   feature$left_sw,
#                                   feature$right_sw,
#                                   subunit.total.intensities,
#                                   subunit.mws,
#                                   subunit.abundances)
#        # Add the resulting information to the initial feature DT.
#        features[i, mw_estimated := res$mw_estimated]
#        features[i, stoichiometry := res$stoichiometry]
#    }
#@CORRECTED looping did result in empty data.table

    # Estimate the stoichiometry of each complex feature
    # as well as its molecular mass.
    res <- lapply(seq(1:nrow(features)), function(i){
       feature=features[i]
       subunits <- strsplit(feature$subgroup, ';')[[1]]
       # Extract only the traces for the subunits that make up this feature.
       subunit.traces <- traces.long[subunits]
       # Extract the relevant information from the protein.info DT.
       subunit.mws <- protein.info[subunits, protein_mw]
       subunit.abundances <- protein.info[subunits, protein_concentration]
       subunit.total.intensities <- protein.info[subunits, total_intensity]
       # Call the helper function to estimate the molecular mass of the
       # complex and its stoichiometry.
       res <- estimateComplexMass(subunit.traces,
         feature$left_sw,
         feature$right_sw,
         subunit.total.intensities,
         subunit.mws,
         subunit.abundances)
    })

    features[, mw_estimated := unlist(lapply(res, function(x) x$mw_estimated))]
    features[, stoichiometry := unlist(lapply(res, function(x) x$stoichiometry))]

    # Produce a score that measures how well both estimates agree.
    features[, mw_delta := abs(mw_estimated - mw_apparent)]

    print(features)  ##### @ERROR THIS IS A BUG IN R I GUESS!!!!!
    return(features)
}






#' A helper function to compute an estimate of the mass of a complex and
#' its stoichiometry.
#' @param traces A long list style data.table holding the intensity
#'     observations of each subunit. This data.table should have the format:
#'     \itemize{
#'      \item \code{protein_id} The subunit identifier
#'      \item \code{fraction} Where the observation was made
#'      \item \code{intensity} The observed intensity
#'     }
#' @param left.boundary The left feature boundary
#' @param right.boundary The right feature boundary
#' @param subunit.total.intensities A numeric vector with the total trace
#'     intensities across all fractions.
#' @param subunit.mws A numeric vector with molecular weights of the subunits.
#' @param subunit.abundances A numeric vector with total abundances of the
#'     subunits.
#' @return A list with the following to entries 'stoichiometry' and
#'     'mw_estimated'.
estimateComplexMass <- function(traces,
                                left.boundary,
                                right.boundary,
                                subunit.total.intensities,
                                subunit.mws,
                                subunit.abundances) {
    subunits <- unique(traces$protein_id)
    n.subunits <- length(subunits)
    stopifnot(length(subunit.mws) == n.subunits &&
              length(subunit.abundances) == n.subunits)
    # For each subunit sum up the intensity of its trace within the
    # feature boundaries.
    # This will produce a named vector of intensities.
    subunit.intensities.within.feature <-
        sapply(subunits, function(subunit) {
            traces[protein_id == subunit &
                   left.boundary <= fraction &
                   fraction <= right.boundary, sum(intensity)]
    })

    # Compute for each subunit how much of its total abundance lies within
    # the boundaries of this feature. This information can then be used
    # to guess the stoichiometry of the complex.
    # Together with the molecular weight of each individual subunit an
    # estimate of the complex mass can be made.
    subunit.relative.intensities <-
        subunit.intensities.within.feature / subunit.total.intensities
    subunit.abundances.in.feature <-
        subunit.relative.intensities * subunit.abundances
    subunit.abundances.in.feature.norm <-
        subunit.abundances.in.feature / min(subunit.abundances.in.feature)
    stoichiometry <- round(subunit.abundances.in.feature.norm)

    est.complex.mass <- sum(stoichiometry * subunit.mws)

    data.table(mw_estimated=est.complex.mass,

               stoichiometry=paste(stoichiometry, collapse=';'))
}
