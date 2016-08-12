#' Reformat result object to a data.table.
#' @param A list containing various results.
#'         \itemize{
#'          \item \code{sw.results} A list of results of the function
#'                \code{findComplexFeatures}. One for each query complex.
#'          \item \code{input.complexes} A character vector of all query
#'                 complexes.
#'          \item \code{corr.cutoff} The correlation cutoff used.
#'          \item \code{window.size} The window size used.
#'          }
#' @return A data.table with all detected complex features of all query complexes
#'          in the format:
#'          \itemize{
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
resultsToTable <- function(swf){
  res <- swf$sw.results
  res.list <- lapply(seq(1:length(res)), function(i){
    features <- res[[i]]$features
    features
  })
  res <- do.call("rbind", res.list)        
  res
}

#' Filter result data.table according to desired complex_ids or minimum correlation score.
#' @param A data.table with the complex features in the format:
#'          \itemize{
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
#'@param complex_ids A character vector containing all desired \code{complex_id} values.
#'@param min_completeness A numeric value specifying the minimum \code{sw_score} for a complex feature.
#'@param collapse_same_peak TRUE or FALSE, collapse if multiple subcomplexes have the same peak apex (identical peak).
#'@return The same data.table format filtered according to the provided parameters.         
subsetComplexFeatures <- function(res,complex_ids=NULL,min_completeness=NULL,collapse_same_peak=FALSE){
  if(!is.null(complex_ids)){
    res <- subset(res,complex_id %in% as.character(complex_ids))
  }
  if(!is.null(min_completeness)){
    res <- subset(res,completeness >= min_completeness)
  }
  if(collapse_same_peak){
    res <- unique(res,by=c("apex",""))
  }
  res
}

#' Reformat result object to a data.table containing only the best subcomplex feature per complex query.
#' @param A list containing various results.
#'         \itemize{
#'          \item \code{sw.results} A list of results of the function
#'                \code{findComplexFeatures}. One for each query complex.
#'          \item \code{input.complexes} A character vector of all query
#'                 complexes.
#'          \item \code{corr.cutoff} The correlation cutoff used.
#'          \item \code{window.size} The window size used.
#'          }
#' @return A data.table with the BEST detected complex feature of each query complexe
#'          in the format:
#'          \itemize{
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
getBestComplexFeature <- function(swf){
  res <- swf$sw.results
  res.list <- lapply(seq(1:length(res)), function(i){
    features <- res[[i]]$features
    features[1] 
  })
  res <- do.call("rbind", res.list)
  res
} 


