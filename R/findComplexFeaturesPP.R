
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
#' @return The same data.table as the input argument extended with the
#'         following columns:
#'         \itemize{
#'          \item \code{n_subunits} The number of subunits in the feature.
#'          \item \code{stoichiometry} The intensity-based stoichiometry.
#'          \item \code{mw_estimated} The estimated molecular weight.
#'          \item \code{mw_apparent} The apparent mw.
#'          \item \code{mw_delta} The delta mw.
#'         }
findComplexFeaturesPP <- function(traces.obj,complexFeaturesSW) {
    features <- complexFeaturesSW$features
    # Compute the number of subunits in each complex feature
    features[, n_subunits := length(strsplit(as.character(subgroup), ';')[[1]]),by=subgroup]

    # create complex trace by summing up intensities
    features.new <- lapply(seq(1:nrow(features)), function(i){
       feature=features[i]
       subunits <- strsplit(feature$subgroup, ';')[[1]]
       # Extract only the traces for the subunits that make up this feature.
       subunit.traces <- subset(traces.obj, subunits)
       subunit.traces.mat <- getIntensityMatrix(subunit.traces)
       # create complex.trace matrix by summing traces for all subunits
       complex.trace <- colSums(subunit.traces.mat)
       # Savitzky-Golay Smoothing (from pracma package)
       complex.trace.SG <-savgol(complex.trace, fl=11, forder = 2, dorder = 0)
       # peak picking (from pracma package)
       complex.peaks <- findpeaks(complex.trace.SG,minpeakdistance=5,nups=3,ndowns=3)
       # convert complex.peaks to data.table
       if (is.null(dim(complex.peaks))){
         complex.peaks <- data.table(intensity=complex.peaks[1],
           apex=complex.peaks[2],
           left_pp=complex.peaks[3],
           right_pp=complex.peaks[4])
       } else {
         complex.peaks <- data.table(complex.peaks)
         setnames(complex.peaks,c("intensity","apex","left_pp","right_pp"))
       }

       # select peaks within boundaries of correlation based window
       ## sel_peaks <- which((complex.peaks$left_pp>=feature$left_sw) & (complex.peaks$right_pp<=feature$right_sw)) # peak boundaries within SW window = problem becaus eof trunkated peaks
       ## sel_peaks <- which((complex.peaks$apex>=feature$left_sw) & (complex.peaks$apex<=feature$right_sw)) # only apex within SW boundaries
       sel_peaks <- which(((complex.peaks$apex>=feature$left_sw) & (complex.peaks$apex<=feature$right_sw)) |
       ((complex.peaks$left_pp>=feature$left_sw) & (complex.peaks$left_pp<feature$right_sw)) |
       ((complex.peaks$right_pp>feature$left_sw) & (complex.peaks$right_pp<=feature$right_sw))) # apex or any peak boundary within SW
       if (length(sel_peaks > 0)) { # peak was detected within SW
         complex.peaks <- complex.peaks[sel_peaks]
         # only peak with highest intensity
         complex.peak <- complex.peaks[which(complex.peaks$intensity==max(complex.peaks$intensity)),]
         complex.peak$area <- sum(complex.trace[complex.peak$left_pp:complex.peak$right_pp])
       } else { #no peak was detected within SW
         #apex=feature$left_sw+((feature$right_sw - feature$left_sw)/2)
         #apex=unique(c(floor(apex),ceiling(apex)))
         #intensity=mean(complex.trace.mat[apex])
         #complex.peak <- data.table(intensity=intensity,apex=mean(apex),left_pp=feature$left_sw,right_pp=feature$right_sw)
         #complex.peak$area <- sum(complex.trace.mat[1,complex.peak$left_pp:complex.peak$right_pp])
         complex.peak <- data.table(intensity=NA,apex=NA,left_pp=NA,right_pp=NA,area=NA)
       }
       complex.peak[,intensity:=NULL]
     })

     features[, apex := unlist(lapply(features.new, function(x) x$apex))]
     features[, left_pp := unlist(lapply(features.new, function(x) x$left_pp))]
     features[, right_pp := unlist(lapply(features.new, function(x) x$right_pp))]
     features[, area := unlist(lapply(features.new, function(x) x$area))]
     data.table(features)

     res <- list(features=features)
     class(res) = 'complexFeaturesPP'

     res
}
