#' Detect subgroups within a window. This is a helper function and should 
#' be called by `findComplexFeatures`.
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
            group.assignment <- c(1, 1)
        } else {
            group.assignment <- c(1, 2)
        }
        return(group.assignment)
    }
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


#' Detect subgroups of proteins within a matrix of protein intensity traces by
#' sliding a window across the SEC dimension. Within each window proteins
#' with traces that correlate well are clustered together.
#'
#' @param trace.mat A numeric matrix where rows correspond to the different
#'        traces.
#' @param protein.names A vector with protein identifiers. This vector has to
#'        have the same length as the number of rows in `trace.mat`.
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
#' # sw.res <- findComplexFeatures(traces, protein.ids)
#' @export
findComplexFeatures <- function(trace.mat,
                                protein.names,
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
    n.zero.entries <- sum(trace.mat == 0)
    noise.mean <- quantile(measure.vals, noise.quantile)
    noise.sd <- sd(measure.vals[measure.vals < noise.mean])
    trace.mat[trace.mat == 0] <- abs(rnorm(n.zero.entries, mean=noise.mean,
                                     sd=noise.sd))
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

    if (nrow(groups.dt) > 0) {
        groups.only <- subset(groups.dt, select=-protein_id)
        setkey(groups.only)
        groups.only <- unique(groups.only)
        groups.only$is_present <- T
        groups.only.wide <-
            cast(groups.only, subgroup ~ sec, value='is_present', fill=F)
        groups.mat <- as.matrix(subset(groups.only.wide, select=-subgroup))
        groups.feats <- findFeatureBoundaries(groups.mat,
                                              groups.only.wide$subgroup,
                                              groups.dt)
        groups.feats <- extendComplexFeatures(groups.feats)
    } else {
        groups.feats <- data.frame()
    }

    result <- list(features=groups.feats,
                   window.size=window.size,
                   corr.cutoff=corr.cutoff)
    class(result) <- 'complexFeaturesSW'
    result
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
#' @return The same data.table as the input argument extended with the
#'         following columns:
#'         \itemize{
#'          \item \code{n_subunits} The number of subunits in the feature.
#'          \item \code{stoichiometry} The intensity-based stoichiometry.
#'          \item \code{mw_estimated} The estimated molecular weight.
#'          \item \code{mw_apparent} The apparent mw.
#'          \item \code{mw_delta} The delta mw.
#'         }
extendComplexFeatures <- function(features) {
    features[, n_subunits := length(strsplit(as.character(subgroup), ';')[[1]]),
             by=subgroup]
    features[, mw_apparent := convertSECToMW(right_sec - left_sec) / 2,
             by=subgroup]
    # TODO
    features[, mw_estimated := 0, by=subgroup]
    features[, mw_delta := 0, by=subgroup]
    features[, stoichiometry := '1A:2B', by=subgroup]
    features
}
