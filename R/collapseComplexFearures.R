#' Collapse detected complex fetures by common apex and subunits.
#'
#' @param complexFeature Object of class ComplexFeaturePP
#' @param rt_height Height at which to cur apex clustering tree.
#' @return An object of type \code{collapsedComplexFeatures} that is a list
#'        containing the following:
#'        \itemize{
#'          \item \code{feature} data.table containing complex feature candidates in the following format:
#'           \itemize{
#'           \item \code{subgroup} The protein_ids of the feature separated by semi-colons.
#'           \item \code{left_sw} The left boundary of the sliding-window feature.
#'           \item \code{right_sw} The right boundary of the sliding-window feature.
#'           \item \code{score} The intra-sliding-window-feature correlation.
#'           }
#'        }

collapseComplexFeatures <- function(complexFeature=complexFeaturesPP,rt_height=5,collapse_method="apex_only"){
  complexFeature <- subset(complexFeature,apex != "NA")
  if (nrow(complexFeature) > 1) {
    apex_dist <- dist(complexFeature$apex)
    apex_clust <- hclust(apex_dist)
    #pdf(paste0("cluster_apex_",gsub("/","",complexFeature$complex_name),".pdf"))
    #  plot(apex_clust)
    #dev.off()
    apex_groups <- cutree(apex_clust, h=rt_height)
    unique_apex_groups <- unique(apex_groups)
    for (apex_group in unique_apex_groups) {
      data <- complexFeature[which(apex_groups == apex_group)]
      if (nrow(data) > 1) {
        complex_subunits <- strsplit(data$subgroup, ';')
        n_complexes <- length(complex_subunits)
        if (collapse_method=="apex_only"){
          #apex_only collapses by apex
          c_proteins <- unique(unlist(complex_subunits))
          c_groups <- rep(1,length(c_proteins))
        } else if (collapse_method == "apex_network"){
          #apex_network collapses by apex and connected network cluster
          bi_list <- lapply(complex_subunits,combn,m=2)
          bi_list.t <- lapply(bi_list,t)
          bi_list.tt <- lapply(bi_list.t,as.data.table)
          binary_interactions <- rbindlist(bi_list.tt)
          binary_interactions <- unique(binary_interactions)
          g <- graph.data.frame(binary_interactions)
          c <- clusters(g)
          c_proteins <- names(c$members)
          c_groups <- as.vector(c$membership)
        }
        unique_c_groups <- unique(c_groups)
        for (c_group in unique_c_groups) {
          tf <- function(a){all(a %in% c_proteins[which(c_groups == c_group)])}
          c_group_member_features <- lapply(complex_subunits,tf)
          c_data <- data[which(c_group_member_features==TRUE)]
          #c_data[,n_primary_features := nrow(c_data)]
          unique_subunits <- unique(unlist(complex_subunits[which(c_group_member_features==TRUE)]))
          c_data[,n_old := n_subunits]
          c_data$n_subunits = length(unique_subunits)
          c_data$subgroup = paste0(unique_subunits,collapse=";")
          c_data <- c_data[order(-n_old,-score,-area)]
          c_data[,n_old := NULL]
          if (exists("new_complexFeature")) {
            new_complexFeature <- rbind(new_complexFeature,c_data[1])
          } else {
            new_complexFeature <- c_data[1]
          }
          rm(c_data)
        }
      } else {
        data <- data[order(-n_subunits,-score,-area)]
        if (exists("new_complexFeature")) {
          new_complexFeature <- rbind(new_complexFeature,data[1])
        } else {
          new_complexFeature <- data[1]
        }
      }
      rm(data)
    }
  } else {
    new_complexFeature <- complexFeature
  }
  return(new_complexFeature[])
}
