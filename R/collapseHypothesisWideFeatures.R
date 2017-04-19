
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

#' getDistance
#' @description getDistance.
#' @param A chracter vector
#' @param B character vector
#' @return numeric distance of subunit overlap
#' @export
getDistance <- function(A,B){
  all_A <- length(A)
  all_B <- length(B)
  unique_A <- length(which(! A %in% B))
  unique_B <- length(which(! B %in% A))
  scoreA <- unique_A/all_A
  scoreB <- unique_B/all_B
  if((scoreA==0) & (scoreB==0)){
    dist <- 0
  } else {
    dist <- 1+(scoreA*scoreB)
  }
  dist
}

#' getDistanceMatrix
#' @description getDistanceMatrix.
#' @param complexFeature data.table containing filtered complex feature results.
#' @return distance matrix
#' @export
getDistanceMatrix <- function(complexFeatures){
  complex_subunits <- strsplit(complexFeatures$subunits_detected, ';')
  n_complexes <- length(complex_subunits)
  dist_mat <- matrix(, nrow = n_complexes, ncol = n_complexes)
  for(m in 1:n_complexes){
    for (n in 1:n_complexes){
      dist_mat[m,n] <- getDistance(complex_subunits[[m]],complex_subunits[[n]])
    }
  }
  dist_mat <- as.dist(dist_mat)
  dist_mat
}

#' complexClustering
#' @description complexClustering.
#' @import ggdendro
#' @param complexFeature data.table containing filtered complex feature results.
#' @param dist_mat distance matrix from getDistanceMatrix
#' @return cluster object
#' @export
complexClustering <- function(complexFeatures,dist_mat){
  hc <- hclust(dist_mat)
  hc$labels = complexFeatures$consecutive_feature_identifier
  hc
}
