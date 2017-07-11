#' getDistance
#' @description Calculate the distance between 2 complex hypotheses.
#' @param A chracter vector
#' @param B character vector
#' @param all_A Numeric, length of chracter vector A.
#' @param all_B Numeric, length of chracter vector B.
#' @return Numeric, distance of subunit overlap.
#' @export
getDistance <- function(A,B, all_A, all_B){
  # all_A <- length(A)
  # all_B <- length(B)
  unique_A <- length(setdiff(A, B))
  unique_B <- length(setdiff(B, A))
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
#' @description calculate the pairwise distance between complex hypotheses.
#' @param complexFeatures data.table containing filtered complex feature results.
#' Must contain the following columns:
#' \itemize{
#' \item subunits_detected: Character, protein ids of complex subunits separated by ';'
#' }
#' @return Numeric distance matrix. 0 = identical; 1 = one complex is a subset of the other;
#' 2 or less = fractional overlap
#' @export
#' @examples
#' ## load example data
#' complexHypotheses <- exampleComplexHypotheses
#' complexHypotheses <- complexHypotheses[, .(subunits_detected = paste0(protein_id, collapse = ";")),
#'                                        by = complex_id] 
#' ## Get the distance matrix
#' distMat <- getDistanceMatrix(complexHypotheses)
#' 
getDistanceMatrix <- function(complexFeatures){
  complex_subunits <- strsplit(complexFeatures$subunits_detected, ';')
  n_complexes <- length(complex_subunits)
  n_subunits <- sapply(complex_subunits,length)
  dist_mat <- matrix(nrow = n_complexes, ncol = n_complexes)
  for(m in 1:n_complexes){
    for (n in 1:n_complexes){
      dist_mat[m,n] <- getDistance(complex_subunits[[m]],complex_subunits[[n]],
                                   n_subunits[m], n_subunits[n])
    }
  }
  dist_mat <- as.dist(dist_mat)
  dist_mat
}

#' complexClustering
#' @description complexClustering.
#' @param complexFeature data.table containing filtered complex feature results.
#' @param dist_mat distance matrix from getDistanceMatrix
#' @return cluster object
#' @export
complexClustering <- function(complexFeatures,dist_mat){
  hc <- hclust(dist_mat)
  hc$labels = complexFeatures$consecutive_feature_identifier
  hc
}


#' Collapse redundant hypotheses
#' @description Remove redundancy in existing complex hypotheses, internal function
#' @import data.table
#' @param hypothesis data.table with complex hypotheses
#' @param redundancy_cutoff numeric maximum overlap distance between two hypotheses (0=identical,1=subset,between 1 and 2=some shared subunits, 2=no shared subunits), default=1
#' @return data.table in the format of complex hypotheses

.collapseWideHypothesis <- function(hypothesis,
                               redundancy_cutoff=1){
  
  dist_hyp <- getDistanceMatrix(hypothesis)
  clust_hyp <- complexClustering(hypothesis,dist_hyp)
  tree_cut=cutree(clust_hyp,h=redundancy_cutoff)
  hypothesis[,consecutive_feature_identifier := .I]
  hypothesis[,unique_feature_identifier := 0]
  unique_feature_identifier=0
  tree_features <- names(tree_cut)
  unique_tree_groups <- unique(tree_cut)
  for (tree_group in unique_tree_groups) {
    unique_feature_identifier = unique_feature_identifier+1
    tree_complex_ids <- which(tree_cut==tree_group)
    hypothesis$unique_feature_identifier[tree_complex_ids] = unique_feature_identifier
  }
  hypothesis <- hypothesis[order(unique_feature_identifier,-n_subunits)]
  if(any(names(hypothesis) == "complex_name")){
      hypothesis[,complex_name := paste(complex_name,collapse="; "),by="unique_feature_identifier"]
  }else{
    hypothesis[,complex_name := paste(complex_id,collapse="_"),by="unique_feature_identifier"]
  }
  hypothesis_unique <- unique(hypothesis,by="unique_feature_identifier")
  hypothesis_unique <- subset(hypothesis_unique,select=c("complex_id","complex_name","subunits_detected"))
  hypothesis_unique <- hypothesis_unique[,list(protein_id = unlist(strsplit(subunits_detected, ";"))), by=c("complex_id","complex_name")]
  return(hypothesis_unique)
}

