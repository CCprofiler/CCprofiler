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
estimateComplexFeatureStoichiometry <- function(traces.obj,complexFeaturesPP) {

    features <- complexFeaturesPP$features

    sel_na <- which(is.na(features$apex))
    if (length(sel_na) > 0) {
      if (length(sel_na) < nrow(features)) {
        features <- features[-sel_na,]
      } else {
        features <- data.frame()
      }
    }


     stoichiometry <- lapply(seq(1:nrow(features)), function(i){
        feature=features[i]
        subunits <- strsplit(feature$subgroup, ';')[[1]]
        fractions <- feature$left_pp:feature$right_pp
        traces_sub <- subset(traces.obj,trace_ids=subunits,fraction_ids=fractions)
        traces_sub.long <- toLongFormat(traces_sub$traces)
        protein.info <- traces_sub.long[, list(total_intensity=sum(intensity)), by=id]
        protein.info[,intensity_ratio:=total_intensity/min(protein.info$total_intensity)]
        protein.info[,stoichiometry:=round(intensity_ratio)]
        data.table(id=paste(protein.info$id, collapse=';'),
        total_intensity=paste(protein.info$total_intensity, collapse=';'),
        intensity_ratio=paste(protein.info$intensity_ratio, collapse=';'),
        stoichiometry=paste(protein.info$stoichiometry, collapse=';'))
      }
      )

      stoichiometry <- do.call("rbind", stoichiometry)
      features <- cbind(features,stoichiometry)
      features[,subgroup:=NULL]

      data.table(features)

      res <- list(features=features)
      class(res) = 'complexFeatureStoichiometry'

      res
   }
