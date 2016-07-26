#.datatable.aware=TRUE

#' Detect subgroups of proteins within a matrix of protein intensity traces by
#' sliding a window across the SEC dimension. Within each window proteins
#' with traces that correlate well are clustered together.
#'
#' @param trace.mat A numeric matrix where rows correspond to the different
#'        traces.
#' @param protein.names A vector with protein identifiers. This vector has to
#'        have the same length as the number of rows in `trace.mat`.
#' @param protein.mw.conc A data.table that stores the molecular weight and
#'        estimate of the absolute abundance for each subunit.
#'        \itemize{
#'         \item \code{protein_id}
#'         \item \code{protein_mw}
#'         \item \code{protein_concentration}
#'        }
#' @param corr.cutoff The correlation value for chromatograms above which
#'        proteins are considered to be coeluting.
#' @param window.size Size of the window. Numeric.
#' @param noise.quantile The quantile to use in estimating the noise level.
#'        Intensity values that are zero are imputed with random noise
#'        according to the noise estimation.
#' @param min.sec The lowest SEC number in the sample.
#' @return An instance of class `complexFeaturesSW`. This is a list with the
#'         following entries:
#'         \itemize{
#'          \item \code{features} A datatable of features. Each feature
#'              can span several SEC fractions.
#'          \item \code{window.size} The window.size used when running this
#'              function.
#'          \item \code{corr.cutoff} The corr.cutoff used when running this
#'              function
#'         }
#' @examples
#' # NOT RUN:
#' # protein.ids <- corum.complex.protein.assoc[complex_id == 181, protein_id]
#' # traces <- subset(protein.traces[protein_id %in% protein.ids],
#' #                  select=-protein_id)
#' # sw.res <- findComplexFeaturesSW(traces, protein.ids, protein.mw.conc)
#' @export
findComplexFeaturesSW <- function(trace.mat,
                                  protein.names,
                                  protein.mw.conc,
                                  corr.cutoff=0.95,
                                  window.size=15,
                                  with.plot=F,
                                  noise.quantile=0.2,
                                  min.sec=1) {
    if (!any(class(trace.mat) == 'matrix')) {
        trace.mat <- as.matrix(trace.mat)
    }
    if (nrow(trace.mat) < 2) {
        stop('Trace matrix needs at least 2 proteins')
    }
    # Impute noise for missing intensity measurements
    measure.vals <- trace.mat[trace.mat != 0]
    #@NEW because of noise imputation error
    if (length(measure.vals) < 10) {
      stop('Trace matrix too sparse for noise imputation')
    }
    n.zero.entries <- sum(trace.mat == 0) # number of ZERO values in matrix
    noise.mean <- quantile(measure.vals, noise.quantile) # mean intensity of noise.quantile (e.g. noise.quantile = lowest 20% of intensities)
    noise.sd <- sd(measure.vals[measure.vals < noise.mean]) # SD of the intensities in the noise.quantile
    trace.mat[trace.mat == 0] <- abs(rnorm(n.zero.entries, mean=noise.mean,
                                     sd=noise.sd)) # replace ZEROs in trace.mat by imputed noise bas don normal distribution
    # Where to stop the sliding window
    end.idx <- ncol(trace.mat) - window.size

    # Analyze each window for groups of protein chromatograms that correlate
    # well.
    groups.by.window <- lapply(seq(1, ncol(trace.mat)), function(i) {
        start.window.idx <- min(end.idx, i)
        end.window.idx <- start.window.idx + window.size
        window.trace.mat <- trace.mat[, start.window.idx:end.window.idx]
        groups.within.window <-
            findComplexFeaturesWithinWindow(window.trace.mat,
                                            corr.cutoff=corr.cutoff,
                                            protein.names=protein.names,
                                            with.plot=with.plot,
                                            sec=start.window.idx)

        groups.within.window
    })

    groups.dt.list <- lapply(1:length(groups.by.window), function(i) {
        # A list of clusters that were found at position `i`.
        subgroups <- groups.by.window[[i]]
        # Produce a list of data.tables, each DT describes a subgroup in long list
        # format.
        subgroups.dt.list <- lapply(subgroups, function(grp) {
            subunits <- grp$subunits
            mean.corr <- grp$mean.corr
            rt.dt <- data.table(sec=i, protein_id=subunits,
                                n_subunits=length(grp),
                                subgroup=paste(subunits, collapse=';'),
                                groupscore=mean.corr)
            # We want to report a feature RT for each position _within_ the
            # window, where the correlation was high enough. So for example
            # if the window at RT == 20 found some subgroup, then the
            # subgroup should be reported for the interval
            # [20, 20 + window.size].
            # To achieve this, we replicate the data.table and each time
            # change the RT value.
            do.call(rbind, lapply(seq(i, min(i + window.size, ncol(trace.mat))),
               function(t) {
                   new.rt.dt <- rt.dt
                   new.rt.dt$sec <- t + (min.sec - 1)
                   new.rt.dt
               })
            )
        })
        do.call(rbind, subgroups.dt.list)
    })
    groups.dt <- do.call(rbind, groups.dt.list)

    if (!is.null(groups.dt)) {
        groups.only <- subset(groups.dt, select=-protein_id)
        setkey(groups.only) # sorts all columns in ascending order
        groups.only <- unique(groups.only)

        groups.only$is_present <- T
        groups.only.wide <- cast(groups.only, subgroup ~ sec, length, value='is_present', fill=F)

        sec.range=ncol(trace.mat)
        subgroup.range=length(unique(groups.only$subgroup))
        groups.only.wide.tf = matrix(FALSE,ncol=sec.range,nrow=subgroup.range)
        groups.only.wide.m = as.matrix(subset(groups.only.wide, select=-subgroup))
        for(n in seq(1,nrow(groups.only.wide.m),1)){
          for(m in seq(1,ncol(groups.only.wide.m),1)){
            if (groups.only.wide.m[n,m] > 0) {  #### @TODO add parameter to function standard 0??
              col=as.numeric(names(subset(groups.only.wide, select=-subgroup))[m])
              groups.only.wide.tf[n,col] = TRUE
            }
          }
        }

        groups.mat <- groups.only.wide.tf
        #groups.mat <- as.matrix(subset(groups.only.wide, select=-subgroup))
        groups.feats <- findFeatureBoundaries(groups.mat,
                                              groups.only.wide$subgroup,
                                              groups.dt)
        protein.mw.conc <- as.data.table(protein.mw.conc)
        #groups.feats <- extendComplexFeatures(groups.feats, trace.mat,
        #                                      protein.names,
        #                                      protein.mw.conc)
        groups.feats <- findFeaturePeaks(groups.feats, trace.mat,
                                              protein.names,protein.mw.conc)
        sel_na <- which(is.na(groups.feats$intensity))
        if (length(sel_na) > 0) {
          if (length(sel_na) < nrow(groups.feats)) {
            groups.feats <- groups.feats[-sel_na,]
          } else {
            groups.feats <- data.frame()
          }
        }
    } else {
        groups.feats <- data.frame()
    }
    # TODO: Include traces object as character string
    # and then get object form environment when plotting.
    result <- list(features=groups.feats,
                   window.size=window.size,
                   corr.cutoff=corr.cutoff)
    class(result) <- 'complexFeaturesSW'
    result
}


#' Detect subgroups within a window. This is a helper function and should
#' be called by `findComplexFeaturesSW`.
#'
#' @param tracemat A matrix of intensity values where chromatograms of
#'        individual proteins are given by the rows.
#' @param corr.cutoff The correlation value for chromatograms above which
#'        proteins are considered to be coeluting.
#' @param protein.names A vector with protein identifiers. This vector has to
#'        have the same length as the number of rows in `tracemat`.
#' @return A list containing lists that describe the found subclusters.
#'         An individual list has the structure:
#'         \itemize{
#'          \item \code{subunits}: character vector
#'          \item \code{mean.corr}: intra cluster correlation
#'         }
findComplexFeaturesWithinWindow <- function(tracemat, corr.cutoff,
                                            protein.names,
                                            with.plot=F, sec=NULL) {
    # In case there are only two proteins in the matrix no clustering
    # is performed.
    if (nrow(tracemat) == 2) {
        corr <- proxy::simil(tracemat, method='correlation')[1]
        if (corr > corr.cutoff) {
            group.assignments <- c(1, 1)
        } else {
            group.assignments <- c(1, 2)
        }
#        return(group.assignment) #@CORRECTED has to be returned as list!
    } else {
    # Compute distance between chromatograms as measured by the pearson
    # correlation.
    distance <- proxy::dist(tracemat, method='correlation')
    # Cluster correlation vectors hierarchically s.t. proteins that correlate
    # well with a similar group of other proteins cluster together.
    cl <- hclust(distance)
    if (with.plot) {
        plot(cl)
        abline(h=1 - corr.cutoff, col='red')
    }
    # Cut the dendrogram at specified distance.
    # For example, if the requested correlation cutoff `corr.cutoff` is
    # 0.7 then the height where the tree is cut is at 1-0.7=0.3.
    # This process will result in a vector of group labels.
    group.assignments <- cutree(cl, h=1 - corr.cutoff)
  }
    # Extract all subclusters and only keep those that have at least 2
    # members. For each subcluster recompute the correlation (without
    # score penalty from imputed noise).
    n.proteins <- nrow(tracemat)
    subclusters.all <- split(1:n.proteins, group.assignments)
    subclusters <- Filter(function(clu) length(clu) >= 2, subclusters.all)
    lapply(subclusters, function(clu) {
        tracemat.clu <- tracemat[clu, ]
        corr.clu <- cor(t(tracemat.clu))
        corr.clu <- corr.clu[upper.tri(corr.clu)]
        corr.clu.mean <- mean(corr.clu)
        list(subunits=protein.names[clu],
             mean.corr=corr.clu.mean
        )
    })
}



# TODO: Implement function to detected feature boundaries using RCPP

#' Helper function to find boundaries of complex features.
#' @param m Logical matrix where an entry m[i, j] indicates if subgroup i was
#' detected at SEC fraction j.
findFeatureBoundaries <- function(m, subgroup.names, groups.dt) {
    boundaries <- lapply(1:nrow(m), function(i) {
        borders <- integer(length=0)
        in.feature <- FALSE
        for (j in 1:ncol(m)) {
            if (m[i, j] && !in.feature && j != ncol(m)) {
                if (m[i, j + 1]) {
                    borders <- c(borders, j)
                    in.feature <- TRUE
                }
            }
            if (!m[i, j] && in.feature) {
                borders <- c(borders, j - 1)
                in.feature <- FALSE
            }
            if (m[i, j] && in.feature && j == ncol(m)) {
                borders <- c(borders, j)
            }
        }
        if (length(borders) > 0) {
            boundaries.m <- matrix(borders, nrow=2)
            left.boundaries <- boundaries.m[1, ]
            right.boundaries <- boundaries.m[2, ]
            subgroup.name <- subgroup.names[i]

            # Extract the groupscore for each feature found
            # for this subgroup.
            n.features.found <- length(left.boundaries)
            do.call(rbind, lapply(1:n.features.found, function(k) {
                right.boundary <- right.boundaries[k]
                left.boundary <- left.boundaries[k]
                groupscore <- groups.dt[subgroup == subgroup.name &
                                        left.boundary <= sec &
                                        sec <= right.boundary,
                                        groupscore][1]
                data.table(subgroup=subgroup.name,
                           left_sec=left.boundaries[k],
                           right_sec=right.boundaries[k],
                           score=groupscore)
            }))
        } else {
            data.table()
        }
    })
    do.call(rbind, boundaries)
}


#' A helper function to extend a list of complex features with additional
#' information.
#'
#' @param features A data.table of complex feature candidates with the
#'        following format:
#'        \itemize{
#'         \item \code{subgroup} A semicolon-separated list of protein
#'                               identifiers.
#'         \item \code{left_sec} The left boundary of the feature.
#'         \item \code{right_sec} The right boundary of the feature.
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
    features[, mw_apparent := convertSECToMW(left_sec+((right_sec - left_sec)/2)),
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
#                                   feature$left_sec,
#                                   feature$right_sec,
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
         feature$left_sec,
         feature$right_sec,
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
