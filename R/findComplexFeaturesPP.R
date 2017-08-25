#' A helper function to perform peak picking on the complex features detected by
#' the sliding-window approach.
#'
#' @param traces.obj An object of type \code{traces.obj}.
#' @param complexFeaturesSW An object of type \code{complexFeaturesSW} that is a list
#'        containing the following:
#'        \itemize{
#'          \item \code{feature} data.table containing complex feature candidates in the following format:
#'           \itemize{
#'           \item \code{subgroup} The protein_ids of the feature separated by semi-colons.
#'           \item \code{left_sw} The left boundary of the sliding-window feature.
#'           \item \code{right_sw} The right boundary of the sliding-window feature.
#'           \item \code{score} The intra-sliding-window-feature correlation.
#'           }
#'          \item \code{window.size} The window.size used when running this
#'              function.
#'          \item \code{corr.cutoff} The corr.cutoff used when running this
#'              function.
#'        }
#' @param smoothing_length Integer length of savitzgy-golay filter.
#' @return An object of type \code{complexFeaturesPP} that is a list
#'        containing the following:
#'        \itemize{
#'          \item \code{feature} data.table containing complex feature candidates in the following format:
#'           \itemize{
#'           \item \code{subgroup} The protein_ids of the feature separated by semi-colons.
#'           \item \code{left_sw} The left boundary of the sliding-window feature.
#'           \item \code{right_sw} The right boundary of the sliding-window feature.
#'           \item \code{score} The intra-sliding-window-feature correlation.
#'           \item \code{n_subunits} The number of protein_ids in the sliding-window feature.
#'           \item \code{apax} The apex of the selected peak by the peak-picker.
#'           \item \code{left_pp} The left boundary of the selected peak by the peak-picker.
#'           \item \code{right_pp} The right boundary of the selected peak by the peak-picker.
#'           \item \code{area} The area (entire complex) of the selected peak by the peak-picker.
#'           }
#'        }

findComplexFeaturesPP <- function(traces.obj,complexFeaturesSW,smoothing_length=11,rt_height=5) {
    features <- copy(complexFeaturesSW)
    # Compute the number of subunits in each complex feature
    features[, n_subunits := length(strsplit(as.character(subgroup), ';')[[1]]),by=subgroup]

    # create complex trace by summing up intensities
    features.new <- lapply(seq(1:nrow(features)), function(i){
      #print(i) for debugging
       feature=features[i]
       subunits <- strsplit(feature$subgroup, ';')[[1]]
       # Extract only the traces for the subunits that make up this feature.
       subunit.traces <- subset(traces.obj, subunits)
       subunit.traces.mat <- getIntensityMatrix(subunit.traces)
       # create complex.trace matrix by summing traces for all subunits
       complex.trace <- colSums(subunit.traces.mat)
       # Savitzky-Golay Smoothing (from pracma package)
       complex.trace.SG <-savgol(complex.trace, fl=smoothing_length, forder = 2, dorder = 0)
       # peak picking (from pracma package)
       complex.peaks <- findpeaks(complex.trace.SG,minpeakdistance=rt_height,nups=3,ndowns=3)
       # convert complex.peaks to data.table
       if (is.null(dim(complex.peaks))){
         if(!is.null(complex.peaks)){
           complex.peaks <- data.table(intensity=complex.peaks[1],
             apex=complex.peaks[2],
             left_pp=complex.peaks[3],
             right_pp=complex.peaks[4])
         } else {
           complex.peaks <- data.table(intensity=NA,
                                       apex=NA,
                                       left_pp=NA,
                                       right_pp=NA)
         }
       } else {
         complex.peaks <- data.table(complex.peaks)
         setnames(complex.peaks,c("intensity","apex","left_pp","right_pp"))
       }

       # select peaks within boundaries of correlation based window
       ## sel_peaks <- which((complex.peaks$left_pp>=feature$left_sw) & (complex.peaks$right_pp<=feature$right_sw)) # peak boundaries within SW window = problem becaus eof trunkated peaks
        sel_peaks <- which((complex.peaks$apex>=feature$left_sw) & (complex.peaks$apex<=feature$right_sw)) # only apex within SW boundaries
       ##sel_peaks <- which(((complex.peaks$apex>=feature$left_sw) & (complex.peaks$apex<=feature$right_sw)) |
       ##((complex.peaks$left_pp>=feature$left_sw) & (complex.peaks$left_pp<feature$right_sw)) |
       ##((complex.peaks$right_pp>feature$left_sw) & (complex.peaks$right_pp<=feature$right_sw))) # apex or any peak boundary within SW
       if (length(sel_peaks > 0)) { # peak was detected within SW
         complex.peaks <- complex.peaks[sel_peaks]
         # only peak with highest intensity
         #complex.peak <- complex.peaks[which(complex.peaks$intensity==max(complex.peaks$intensity)),]
         #complex.peak$area <- sum(complex.trace[complex.peak$left_pp:complex.peak$right_pp])
         complex.peak <- complex.peaks[,area:=NA]
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

      features_rep <- unlist(lapply(features.new, nrow))
      features <- features[rep(seq(1,nrow(features),1),features_rep)]
      features[, apex := unlist(lapply(features.new, function(x) x$apex))]
      features[, left_pp := unlist(lapply(features.new, function(x) x$left_pp))]
      features[, right_pp := unlist(lapply(features.new, function(x) x$right_pp))]
      features[, area := unlist(lapply(features.new, function(x) x$area))]
      data.table(features)

      res <- list(features=features)
      class(res) = 'complexFeaturesPP'

     res
}
