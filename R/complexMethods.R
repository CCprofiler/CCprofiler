#' Filter complex feature table
#' @description Filter result data.table according to desired complex_ids or minimum correlation score.
#' @param res A data.table with the complex features.
#' @param complex_ids A character vector containing all desired \code{complex_id} values.
#' @param min_completeness Numeric between 0 and 1, specifying the required completeness of a feature reltive to the tested hypothesis (keeps all features if at least one is bigger than the cutoff).
#' @param min_subunits Integer specifying minimum number of subunits in a complex.
#' @param min_peak_corr Numeric value betwee 0 and 1 specifying minimum peak correlation, default is 0.5
#' @param min_monomer_distance_factor Numeric value specifying factor to multiply largest monomer MW, default is 1
#' @return The same data.table format filtered according to the provided parameters.
#' @export
subsetComplexFeatures <- function(res,complex_ids=NULL,min_completeness=NULL,min_subunits=NULL,min_peak_corr=0.5,min_monomer_distance_factor=1){
  if(!is.null(complex_ids)){
    res <- subset(res,complex_id %in% as.character(complex_ids))
  }
  if(!is.null(min_completeness)){
    # res <- subset(res,completeness >= min_completeness)
    allowed_ids <- res[completeness>=min_completeness, unique(complex_id)]
    res <- res[complex_id %in% allowed_ids]
  }
  if(!is.null(min_subunits)){
    res <- subset(res,n_subunits_detected >= min_subunits)
  }
  if(!is.null(min_peak_corr)){
    res <- subset(res,peak_corr >= min_peak_corr)
  }
  if(!is.null(min_monomer_distance_factor)){
    dist = lapply(seq(1:nrow(res)), function(i){
      feature=res[i]
      max_monomer_mw <- max(as.numeric(strsplit(feature$monomer_mw, ';')[[1]]))
      mw_dist <- feature$apex_mw-(min_monomer_distance_factor*max_monomer_mw)
      mw_dist
    })
    dist_vector <- unlist(dist)
    sel <- which(dist<=0)
    if(length(sel)>0){
      res <- res[-sel]
    }
  }
  res
}


#' Filter protein feature table
#' @description Filter result data.table according to desired protein_ids or minimum correlation score.
#' @param res A data.table with the protein features.
#' @param protein_ids A character vector containing all desired \code{protein_id} values.
#' @param min_completeness Numeric between 0 and 1, specifying the required completeness of a feature reltive to the tested hypothesis.
#' @param min_subunits Integer specifying minimum number of subunits in a protein.
#' @param min_peak_corr Numeric value betwee 0 and 1 specifying minimum peak correlation, default is 0.5
#' @return The same data.table format filtered according to the provided parameters.
#' @export
subsetProteinFeatures <- function(res,protein_ids=NULL,min_completeness=NULL,min_subunits=NULL,min_peak_corr=0.5){
  if(!is.null(protein_ids)){
    res <- subset(res,protein_id %in% as.character(protein_ids))
  }
  if(!is.null(min_completeness)){
    res <- subset(res,completeness >= min_completeness)
  }
  if(!is.null(min_subunits)){
    res <- subset(res,n_subunits_detected >= min_subunits)
  }
  if(!is.null(min_peak_corr)){
    res <- subset(res,peak_corr >= min_peak_corr)
  }
  res
}
