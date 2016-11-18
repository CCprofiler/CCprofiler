#' A helper function to extend a list of complex features with additional
#' stoichiometry estimation.
#'
#' @param traces.obj An object of type \code{traces.obj}.
#' @param complexFeaturesPP An object of type \code{complexFeaturesPP} that is a list
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
#' @return An object of type \code{complexFeatureStoichiometries} that is a list
#'        containing the following:
#'        \itemize{
#'          \item \code{feature} data.table containing complex feature candidates in the following format:
#'           \itemize{
#'           \item \code{left_sw} The left boundary of the sliding-window feature.
#'           \item \code{right_sw} The right boundary of the sliding-window feature.
#'           \item \code{score} The intra-sliding-window-feature correlation.
#'           \item \code{n_subunits} The number of protein_ids in the sliding-window feature.
#'           \item \code{apax} The apex of the selected peak by the peak-picker.
#'           \item \code{left_pp} The left boundary of the selected peak by the peak-picker.
#'           \item \code{right_pp} The right boundary of the selected peak by the peak-picker.
#'           \item \code{area} The area (entire complex) of the selected peak by the peak-picker.
#'           \item \code{id} The protein_ids of the feature separated by semi-colons.
#'           \item \code{total_intensity} The intensity of all protein_ids of the feature separated by semi-colons.
#'           \item \code{intensity_ratio} The intensity ratio of all protein_ids of the feature separated by semi-colons.
#'           \item \code{stoichiometry} The rounded \code{intensity_ratio} of all protein_ids of the feature separated by semi-colons.
#'           }
#'        }
#' @export

estimateComplexFeatureStoichiometry <- function(traces.obj,complexFeaturesPP) {

  features <- complexFeaturesPP$features
  # remove sliding-window features where no peak was detected (these have an NA in the features$apex)
  sel_na <- which(is.na(features$apex))
  if (length(sel_na) > 0) {
    if (length(sel_na) < nrow(features)) {
      features <- features[-sel_na,]
    } else {
      features <- data.frame()
    }
  }

  # estimate the stoichiometry of each detected feature with a picked peak
  stoichiometry <- lapply(seq(1:nrow(features)), function(i){
    feature=features[i]
    subunits <- strsplit(feature$subgroup, ';')[[1]]
    # select sec fractions whithin the boundaries of the picked peak and subset the traces.obj
    fractions <- feature$left_pp:feature$right_pp
    traces_sub <- subset(traces.obj,trace_ids=subunits,fraction_ids=fractions)
    traces_sub.long <- toLongFormat(traces_sub$traces) #long format is easier for processing
    # make sure intensity is numeric
    traces_sub.long$intensity <- as.numeric(traces_sub.long$intensity)
    # sum up the intensities for each subunit within the peak boundaries
    protein.info <- traces_sub.long[, list(total_intensity=sum(intensity)), by=id]
    # build intensity ratio by dividing all total subunit intensities by the value of the subunit with the lowest intensity
    protein.info[,intensity_ratio:=total_intensity/min(protein.info$total_intensity)]
    # round the intensity ratios to get integer stoichiometry estimates
    protein.info[,stoichiometry:=round(intensity_ratio)]
    data.table(id=paste(protein.info$id, collapse=';'),
                total_intensity=paste(protein.info$total_intensity, collapse=';'),
                intensity_ratio=paste(protein.info$intensity_ratio, collapse=';'),
                stoichiometry=paste(protein.info$stoichiometry, collapse=';'))
  }
  )

  stoichiometry <- do.call("rbind", stoichiometry)
  features <- cbind(features,stoichiometry)
  # remove subgroup column becaus ewe have new id column which has the same protein_id order as the intensities, ratios and stoichiometry estimates
  features[,subgroup:=NULL]

  data.table(features)
  res <- list(features=features)
  class(res) = 'complexFeatureStoichiometry'
  res
}
