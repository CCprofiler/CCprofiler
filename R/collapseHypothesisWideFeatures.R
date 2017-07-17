#' Determine unique feature groups
#' @description Determine unique feature groups based on the apex distances in the fraction dimension 
#' as well as the overlap between co-eluting subunits.
#' @param feature_table data.table containing filtered complex feature results.
#' @param rt_height Numeric rt difference for clustering along the fractionation dimension, default = 2.
#' @param distance_cutoff Numeric distnce of detected subunit overlap (0=identical, 1=subcomplexes, 2=no shared subunits), default = 1.
#' @param distance_cutoff Numeric distnce of detected subunit overlap (0=identical, 1=subcomplexes, 2=no shared subunits), default = 1.
#' @param plot Logical whether to plot the fraction dimenson apex distance histogram, default = \code{FALSE}.
#' @return data.table containing filtered complex feature results with an extra column with a unique_feature_identifier.
#' @export
#' @examples 
#' ## run unique feature annotation
#' complexFeaturesUnique <- getUniqueFeatureGroups(exampleComplexFeatures,
#'                                                 rt_height=2,
#'                                                 distance_cutoff=1,
#'                                                 plot=TRUE)
#'                                                 
getUniqueFeatureGroups <- function(feature_table,
                                   rt_height=2,
                                   distance_cutoff=1,
                                   plot=FALSE){
  apex_dist <- dist(feature_table$apex)
  apex_clust <- hclust(apex_dist)
  if (plot) {
    plot(apex_clust)
  }
  apex_groups <- cutree(apex_clust, h=rt_height)
  unique_apex_groups <- unique(apex_groups)
  feature_table[,consecutive_feature_identifier := .I]
  feature_table[,unique_feature_identifier := 0]
  unique_feature_identifier = 0
  for (apex_group in unique_apex_groups) {
    data <- feature_table[which(apex_groups == apex_group)]
    if (nrow(data) > 1) {
      dist <- getDistanceMatrix(data)
      clust <- complexClustering(complexFeatures=data,dist_mat=dist)
      tree_cut=cutree(clust,h=distance_cutoff)
      tree_features <- names(tree_cut)
      unique_tree_groups <- unique(tree_cut)
      for (tree_group in unique_tree_groups) {
        unique_feature_identifier = unique_feature_identifier+1
        tree_complex_ids <- names(tree_cut)[which(tree_cut==tree_group)]
        feature_table$unique_feature_identifier[which(feature_table$consecutive_feature_identifier %in% tree_complex_ids)] = unique_feature_identifier
      }
    } else {
      unique_feature_identifier = unique_feature_identifier+1
      feature_table$unique_feature_identifier[which(feature_table$consecutive_feature_identifier == data$consecutive_feature_identifier)] = unique_feature_identifier
    }
    rm(data)
  }
  feature_table[,consecutive_feature_identifier := NULL]
  feature_table <- feature_table[order(unique_feature_identifier,-n_subunits_detected,-sw_score,-area,mw_diff)]
  return(feature_table)
}

#' Collapse features by unique feature groups
#' @description Collapse features with the same unique_feature_identifier as determined by \code{\link{getUniqueFeatureGroups}}.
#' @param feature_table data.table containing filtered complex feature results with an extra column with a unique_feature_identifier.
#' @param rm_decoys Logical whether to remove decoys before collapsing, default = \code{TRUE}.
#' @return data.table containing collapsed complex features.
#' @export
#' @examples 
#' ## run unique feature annotation
#' complexFeaturesUnique <- getUniqueFeatureGroups(exampleComplexFeatures,
#'                                                 rt_height=2,
#'                                                 distance_cutoff=1,
#'                                                 plot=TRUE)
#' ## Collapse across unique features
#' complexFeaturesCollapsed <- callapseByUniqueFeatureGroups(complexFeaturesUnique,
#'                                                           rm_decoys=TRUE)
#' ## Inspect output
#' complexFeaturesCollapsed
#' 
callapseByUniqueFeatureGroups <- function(feature_table,rm_decoys=TRUE){
  features <- copy(feature_table)
  if(rm_decoys){
    decoy_ind <- grepl("DECOY",features$complex_id)
    features <- features[!decoy_ind]
  }
  unique_feature_groups <- unique(features$unique_feature_identifier)
  collapsedFeatures <- rbindlist(lapply(unique_feature_groups,.collapseWithinFeatureGroup,feature_table=features))
  return(collapsedFeatures)
}


.collapseWithinFeatureGroup <- function(unique_feature_group, feature_table){
  features <- subset(feature_table,unique_feature_identifier==unique_feature_group)
  collapsed <- data.table(
    complex_id = paste(features$complex_id,collapse=";"),
    complex_name = paste(features$complex_name,collapse=";"),
    subunits = paste(unique(unlist(strsplit(features$subunits_detected, ';'))),collapse=";"),
    left_pp = median(features$left_pp),
    right_pp = median(features$right_pp),
    apex = median(features$apex),
    unique_feature_identifier = unique(features$unique_feature_identifier))
  return(collapsed)
}