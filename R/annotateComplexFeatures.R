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
annotateComplexFeatures <- function(traces.obj,complexFeatureStoichiometries,complex.annotation) {
  
  setkey(complex.annotation, protein_id)
  
  features <- complexFeatureStoichiometries$features
  features[,complex_id := complex.annotation$complex_id[1]]
  features[,complex_name := complex.annotation$complex_name[1]]
  features[,n_subunits_annotated := nrow(complex.annotation)]
  features[,completeness := n_subunits/n_subunits_annotated]
  features[,subunits_annotated := paste(complex.annotation$protein_id, collapse=';')]
  
  protein.mw <- subset(traces.obj$trace_annotation,id %in% complex.annotation$protein_id)
  setkey(protein.mw, id)
  
  mw <- lapply(seq(1:nrow(features)), function(i){
    feature=features[i]
    subunits <- strsplit(feature$id, ';')[[1]]
    
    subunit_MW <-  protein.mw$protein_mw[protein.mw$id %in% subunits]
    
    subunit_SEC <- (log(subunit_MW)-9.682387)/(-0.1043329) # @TODO function
    
    stoichiometry <- strsplit(feature$stoichiometry, ';')[[1]]
    stoichiometry <- as.integer(stoichiometry)
    
    complex_mw <- sum(stoichiometry*subunit_MW)
    
    complex_SEC <-  (log(complex_mw)-9.682387)/(-0.1043329)
    
    SEC_diff <- abs(feature$apex - complex_SEC)
    
    data.table(monomer_mw=paste(subunit_MW, collapse=';'),
               monomer_sec=paste(subunit_SEC, collapse=';'),
               complex_mw_estimated=complex_mw,
               complex_sec_estimated=complex_SEC,
               sec_diff=SEC_diff)
  }
  )
  
  mw <- do.call("rbind", mw)
  features <- cbind(features,mw)
  
  setcolorder(features, c("complex_id", "complex_name", "subunits_annotated",
                   "n_subunits_annotated","id","n_subunits",
                   "completeness","left_sw","right_sw","score",
                   "left_pp","right_pp","apex","area",
                   "total_intensity","intensity_ratio","stoichiometry",
                   "monomer_mw","monomer_sec","complex_mw_estimated","complex_sec_estimated","sec_diff"))
  
  setnames(features,c("complex_id", "complex_name", "subunits_annotated",
                      "n_subunits_annotated","subunits_detected","n_subunits_detected",
                      "completeness","left_sw","right_sw","sw_score",
                      "left_pp","right_pp","apex","area",
                      "total_intensity","intensity_ratio","stoichiometry_estimated",
                      "monomer_mw","monomer_sec","complex_mw_estimated","complex_sec_estimated","sec_diff"))
  
  features <- features[order(-n_subunits_detected, -area, -sw_score)]
  
  data.table(features)
  
  res <- list(features=features)
  class(res) = 'complexFeaturesAnnotated'
  
  res
}