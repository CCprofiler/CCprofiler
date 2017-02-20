#.datatable.aware=TRUE

#' Detect subgroups of proteins within a matrix of protein intensity traces by
#' sliding a window across the SEC dimension. Within each window proteins
#' with traces that correlate well are clustered together.
#'
#' @param trace.mat A matrix of all traces in a complex hypothesis with small perturbations
#'        added to all 0 values for correlation calculation.
#' @param corr.cutoff The correlation value for chromatograms above which
#'        proteins are considered to be coeluting.
#' @param window_size Size of the window. Numeric.
#' @param with.plot T (TRUE) or F (FALSE) whether to plot the correlation tree.
#' @param noise.quantile The quantile to use in estimating the noise level.
#'        Intensity values that are zero are imputed with random noise
#'        according to the noise estimation.
#' @param min.sec The lowest SEC number in the sample.
#' @return An object of type \code{complexFeaturesSW} that is a list
#'        containing the following:
#'        \itemize{
#'          \item \code{feature} data.table containing complex feature candidates in the following format:
#'           \itemize{
#'           \item \code{subgroup} The protein_ids of the feature separated by semi-colons.
#'           \item \code{left_sw} The left boundary of the sliding-window feature.
#'           \item \code{right_sw} The right boundary of the sliding-window feature.
#'           \item \code{score} The intra-sliding-window-feature correlation.
#'           }
#'          \item \code{window_size} The window_size used when running this
#'              function.
#'          \item \code{corr.cutoff} The corr.cutoff used when running this
#'              function.
#'        }
#' @examples
#' # NOT RUN:
#' # protein.ids <- corum.complex.protein.assoc[complex_id == 181, protein_id]
#' # traces <- subset(protein.traces[protein_id %in% protein.ids],
#' #                  select=-protein_id)
#' # sw.res <- findComplexFeaturesSW(traces, protein.ids, protein.mw.conc)
#' @export
findComplexFeaturesSW <- function(trace.mat,
                                  corr.cutoff=0.95,
                                  window_size=15,
                                  with.plot=F,
                                  min.sec=1){
  # trace.mat = getIntensityMatrix(traces.obj)
  # protein.names = traces.obj$traces$id
  protein.names = rownames(trace.mat)
  #protein.mw.conc <- protein.mw.conc
  if (!any(class(trace.mat) == 'matrix')) {
    trace.mat <- as.matrix(trace.mat)
  }
  if (nrow(trace.mat) < 2) {
    stop('Trace matrix needs at least 2 proteins')
  }
  # # Impute noise for missing intensity measurements
  # measure.vals <- trace.mat[trace.mat != 0]
  # if (length(measure.vals) < 10) {
  #   stop('Trace matrix too sparse for noise imputation')
  # }
  # n.zero.entries <- sum(trace.mat == 0) # number of ZERO values in matrix
  # noise.mean <- quantile(measure.vals, noise.quantile) # mean intensity of noise.quantile (e.g. noise.quantile = lowest 20% of intensities)
  # noise.sd <- sd(measure.vals[measure.vals < noise.mean]) # SD of the intensities in the noise.quantile
  # set.seed(123) # set seed to always get same results
  # trace.mat[trace.mat == 0] <- abs(rnorm(n.zero.entries, mean=noise.mean,
  #                                        sd=noise.sd)) # replace ZEROs in trace.mat by imputed noise bas don normal distribution
  # Where to stop the sliding window
  end.idx <- ncol(trace.mat) - window_size
  
  # Analyze each window for groups of protein chromatograms that correlate
  # well.
  groups.by.window <- lapply(seq(1, ncol(trace.mat)), function(i) {
    start.window.idx <- min(end.idx, i)
    end.window.idx <- start.window.idx + window_size
    window.trace.mat <- trace.mat[, start.window.idx:end.window.idx]
    groups.within.window <-
      findComplexFeaturesWithinWindow(window.trace.mat=window.trace.mat,
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
      # [20, 20 + window_size].
      # To achieve this, we replicate the data.table and each time
      # change the RT value.
      do.call(rbind, lapply(seq(i, min(i + window_size, ncol(trace.mat))),
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
  } else {
    groups.feats <- data.frame()
  }
  # TODO: Include traces object as character string
  # and then get object form environment when plotting.
  result <- list(features=groups.feats,
                 window_size=window_size,
                 corr.cutoff=corr.cutoff)
  class(result) <- 'complexFeaturesSW'
  result
}


#' Detect subgroups within a window. This is a helper function and should
#' be called by `findComplexFeaturesSW`.
#'
#' @param window.trace.mat A matrix of intensity values where chromatograms of
#'        individual proteins are given by the rows.
#' @param protein.names A vector with protein identifiers. This vector has to
#'        have the same length as the number of rows in `window.trace.mat`.
#' @param corr.cutoff The correlation value for chromatograms above which
#'        proteins are considered to be coeluting.
#' @return A list containing lists that describe the found subclusters.
#'         An individual list has the structure:
#'         \itemize{
#'          \item \code{subunits}: character vector
#'          \item \code{mean.corr}: intra cluster correlation
#'         }
findComplexFeaturesWithinWindow <- function(window.trace.mat, protein.names,
                                            corr.cutoff,
                                            with.plot=F, sec=NULL) {
  # In case there are only two proteins in the matrix no clustering
  # is performed.
  if (nrow(window.trace.mat) == 2) {
    corr <- proxy::simil(window.trace.mat, method='correlation')[1]
    if (corr > corr.cutoff) {
      group.assignments <- c(1, 1)
    } else {
      group.assignments <- c(1, 2)
    }
    #        return(group.assignment) #@CORRECTED has to be returned as list!
  } else {
    # Compute distance between chromatograms as measured by the pearson
    # correlation.
    distance <- proxy::dist(window.trace.mat, method='correlation')
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
  n.proteins <- nrow(window.trace.mat)
  subclusters.all <- split(1:n.proteins, group.assignments)
  subclusters <- Filter(function(clu) length(clu) >= 2, subclusters.all)
  lapply(subclusters, function(clu) {
    window.trace.mat.clu <- window.trace.mat[clu, ]
    corr.clu <- cor(t(window.trace.mat.clu))
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
#' @param subgroup.names A vector with subgroup identifiers (protein
#'        identifies separated gy ";"). This vector has to
#'        have the same length as the number of rows in `m`.
#' @param groups.dt A data.table with following entries:
#'        \itemize{
#'         \item \code{sec}
#'         \item \code{protein_id}
#'         \item \code{n_subunits}
#'         \item \code{subgoup}
#'         \item \code{groupscore}
#'        }
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
                   left_sw=left.boundaries[k],
                   right_sw=right.boundaries[k],
                   score=groupscore)
      }))
    } else {
      data.table()
    }
  })
  do.call(rbind, boundaries)
}
