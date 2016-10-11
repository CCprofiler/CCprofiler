#' A helper function to extend a list of complex features with additional
#' information.
#'

#' @param traces.obj An object of type \code{traces.obj}.
#' @param complexFeatureStoichiometries An object of type \code{complexFeatureStoichiometries} that is a list
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
#' @param complex.annotation A data.table that stores the complex name and protein subunits for each complex_id.
#'        \itemize{
#'         \item \code{complex_id}
#'         \item \code{complex_name}
#'         \item \code{protein_id}
#'        }
#' @param MWSECcalibrationFunctions A list that stores the functions for converting SEC fractions to molecular weight and vice versa. This
#'        object is generated from the calibrateSECMW function.
#'        \itemize{
#'         \item \code{MWtoSECfraction} Function that needs MW as input and reports according SEC fraction.
#'         \item \code{SECfractionToMW} Function that needs SEC fraction as input and reports according MW.
#'        }
#' @return An object of type \code{complexFeaturesAnnotated} that is a list
#'        containing the following:
#'        \itemize{
#'          \item \code{feature} data.table containing complex feature candidates in the following format:
#'           \itemize{
#'           \item \code{complex_id} The complex_id of the query complex.
#'           \item \code{complex_name} The complex_name of the query complex.
#'           \item \code{subunits_annotated} The subunits (protein_ids) annotated for the complex separated by semi-colons.
#'           \item \code{n_subunits_annotated} The number of subunits (protein_ids) annotated for the complex separated by semi-colons.
#'           \item \code{subunits_with_signal} The subunits (protein_ids) with an MS/MS signal for the complex separated by semi-colons.
#'           \item \code{n_subunits_with_signal} The number of subunits (protein_ids) with an MS/MS signal for the complex separated by semi-colons.
#'           \item \code{subunits_detected} The subunits (protein_ids) detceted in the feature for the complex separated by semi-colons.
#'           \item \code{n_subunits_detected} The number of subunits (protein_ids) detceted in the feature for the complex separated by semi-colons.
#'           \item \code{completeness} The complex completeness as defined by \code{n_subunits_detected} divided by \code{n_subunits_annotated}.
#'           \item \code{left_sw} The left boundary of the sliding-window feature.
#'           \item \code{right_sw} The right boundary of the sliding-window feature.
#'           \item \code{sw_score} The intra-sliding-window-feature correlation.
#'           \item \code{left_pp} The left boundary of the selected peak by the peak-picker.
#'           \item \code{right_pp} The right boundary of the selected peak by the peak-picker.
#'           \item \code{apax} The apex of the selected peak by the peak-picker.
#'           \item \code{area} The area (entire complex) of the selected peak by the peak-picker.
#'           \item \code{total_intensity} The intensity of all protein_ids of the feature separated by semi-colons.
#'           \item \code{intensity_ratio} The intensity ratio of all protein_ids of the feature separated by semi-colons.
#'           \item \code{stoichiometry_estimated} The rounded \code{intensity_ratio} of all protein_ids of the feature separated by semi-colons.
#'           \item \code{monomer_mw} The monomer molecular weights of all protein_ids of the feature separated by semi-colons.
#'           \item \code{monomer_sec} The monomer sec fraction of all protein_ids of the feature separated by semi-colons.
#'           \item \code{complex_mw_estimated} The complex molecular weight as expected fro the \code{stoichiometry_estimated}.
#'           \item \code{complex_sec_estimated} The complex sec fraction as expected fro the \code{stoichiometry_estimated}.
#'           \item \code{sec_diff} Difference between \code{complex_sec_estimated} and \code{apax} of the feature.
#'          }
#'        }
#' @export

annotateComplexFeatures <- function(traces.obj,complexFeatureStoichiometries,complex.annotation,MWSECcalibrationFunctions) {

  setkey(complex.annotation, protein_id)

  # annotate feature data.table with known complex information
  features <- complexFeatureStoichiometries$features
  features[,complex_id := complex.annotation$complex_id[1]]
  features[,complex_name := complex.annotation$complex_name[1]]
  features[,n_subunits_annotated := nrow(complex.annotation)]
  features[,completeness := n_subunits/n_subunits_annotated]
  features[,subunits_annotated := paste(complex.annotation$protein_id, collapse=';')]

  # extract protein molecular weights from the trace_annotations in the traces.obj
  protein.mw <- subset(traces.obj$trace_annotation,id %in% complex.annotation$protein_id)
  setkey(protein.mw, id)

  # add molecular weight and sec fraction information to each feature
  mw <- lapply(seq(1:nrow(features)), function(i){
    feature=features[i]
    # annotate features by complex information and monomer molecular weights
    subunits_annotated <- strsplit(feature$subunits_annotated, ';')[[1]]
    subunits_with_signal <- traces.obj$trace_annotation$id[i=which(traces.obj$trace_annotation$id %in% subunits_annotated)]
    subunits_with_signal <- sort(subunits_with_signal)
    n_subunits_with_signal <- length(subunits_with_signal)
    subunits <- strsplit(feature$id, ';')[[1]]
    subunit_MW <-  protein.mw$protein_mw[protein.mw$id %in% subunits]
    #subunit_SEC <- (log(subunit_MW)-9.682387)/(-0.1043329) 
    subunit_SEC <- MWSECcalibrationFunctions$MWtoSECfraction(subunit_MW)
    # calculate complex molecular weight
    stoichiometry <- strsplit(feature$stoichiometry, ';')[[1]]
    stoichiometry <- as.integer(stoichiometry)
    complex_mw <- sum(stoichiometry*subunit_MW)
    #complex_SEC <-  (log(complex_mw)-9.682387)/(-0.1043329)
    complex_SEC <- MWSECcalibrationFunctions$MWtoSECfraction(complex_mw)
    # calculate apex molecular weight
    #apex_MW <-  exp((-0.1043329 * feature$apex) + 9.682387)
    apex_MW <- MWSECcalibrationFunctions$SECfractionToMW(feature$apex)
    # calculate difference between apex of selected peak and the estimated comples sec fraction
    SEC_diff <- abs(feature$apex - complex_SEC)
    MW_diff <- abs(apex_MW - complex_mw)
    # create output data.table
    data.table(monomer_mw=paste(subunit_MW, collapse=';'),
               monomer_sec=paste(subunit_SEC, collapse=';'),
               complex_mw_estimated=complex_mw,
               complex_sec_estimated=complex_SEC,
               sec_diff=SEC_diff,
               mw_diff=MW_diff,
               subunits_with_signal = paste(subunits_with_signal, collapse=';'),
               n_subunits_with_signal = n_subunits_with_signal,
               apex_mw = apex_MW)
  }
  )

  mw <- do.call("rbind", mw)
  features <- cbind(features,mw)

  # order features data.table in a usefull way
  setcolorder(features, c("complex_id", "complex_name", "subunits_annotated",
                   "n_subunits_annotated","subunits_with_signal","n_subunits_with_signal",
                   "id","n_subunits",
                   "completeness","left_sw","right_sw","score",
                   "left_pp","right_pp","apex","apex_mw","area",
                   "total_intensity","intensity_ratio","stoichiometry",
                   "monomer_mw","monomer_sec","complex_mw_estimated",
                   "complex_sec_estimated","sec_diff","mw_diff"))

  # provide better column names for features data.table
  setnames(features,c("complex_id", "complex_name", "subunits_annotated",
                      "n_subunits_annotated","subunits_with_signal","n_subunits_with_signal",
                      "subunits_detected","n_subunits_detected",
                      "completeness","left_sw","right_sw","sw_score",
                      "left_pp","right_pp","apex","apex_mw","area",
                      "total_intensity","intensity_ratio","stoichiometry_estimated",
                      "monomer_mw","monomer_sec","complex_mw_estimated",
                      "complex_sec_estimated","sec_diff","mw_diff"))

  # sort the features by the number of detected subunits, the sliding-windoe correlation, and the peak area
  features <- features[order(-n_subunits_detected,-sw_score,-area,mw_diff)]

  data.table(features)
  res <- list(features=features)
  class(res) = 'complexFeaturesAnnotated'
  res
}
