
#' collapseHypothesisWideFeatures
#' @description collapseHypothesisWideFeatures.
#' @param complexFeature data.table containing filtered complex feature results.
#' @param rt_height numeric rt difference for clustering along the sec dimension, default = 2
#' @param distance_cutoff numeric distnce of detected subunit overlap (0=identical, 1=subcomplexes, 2=no shared subunits), default=1
#' @return data.table containing collapsed complex feature results
#' @export
collapseHypothesisWideFeatures <- function(complexFeature=complexFeatures,rt_height=2,distance_cutoff=1){
    apex_dist <- dist(complexFeature$apex)
    apex_clust <- hclust(apex_dist)
    #pdf(paste0("cluster_apex_",gsub("/","",complexFeature$complex_name),".pdf"))
    #  plot(apex_clust)
    #dev.off()
    apex_groups <- cutree(apex_clust, h=rt_height)
    unique_apex_groups <- unique(apex_groups)
    complexFeature[,consecutive_feature_identifier := .I]
    complexFeature[,unique_feature_identifier := 0]
    unique_feature_identifier = 0
    for (apex_group in unique_apex_groups) {
      data <- complexFeature[which(apex_groups == apex_group)]
      if (nrow(data) > 1) {
        dist <- getDistanceMatrix(data)
        clust <- complexClustering(complexFeatures=data,dist_mat=dist)
        tree_cut=cutree(clust,h=distance_cutoff)
        tree_features <- names(tree_cut)
        unique_tree_groups <- unique(tree_cut)
        for (tree_group in unique_tree_groups) {
          unique_feature_identifier = unique_feature_identifier+1
          tree_complex_ids <- names(tree_cut)[which(tree_cut==tree_group)]
          complexFeature$unique_feature_identifier[which(complexFeature$consecutive_feature_identifier %in% tree_complex_ids)] = unique_feature_identifier
        }
      } else {
        unique_feature_identifier = unique_feature_identifier+1
        complexFeature$unique_feature_identifier[which(complexFeature$consecutive_feature_identifier == data$consecutive_feature_identifier)] = unique_feature_identifier
      }
      rm(data)
    }
  complexFeature[,consecutive_feature_identifier := NULL]
  complexFeature <- complexFeature[order(unique_feature_identifier,-n_subunits_detected,-sw_score,-area,mw_diff)]
  complexFeature
}

