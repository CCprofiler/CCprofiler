#' getBestParameterStats_constraintFDR
#' @description getBestParameterStats_constraintFDR. Do x-fold cross-FDR estimation.
#' getBestParameterStats_constraintFDR(complex_feature_grid_filtered=var_complex_feature_grid_filtered,complex_hypotheses=var_complex_hypotheses,n_subsets=var_n_subsets,FDR_cutoff=var_FDR_cutoff)
#' @param complex_feature_grid_filtered list with all complex features from grid search
#' @param complex_hypotheses data.table with complex hypotheses including decoys
#' @param n_subsets numeric number of subsets used for cross-calculation of FDR, default = 3
#' @param FDR_cutoff numeric maximum FDR that should be considered, default = 0.1
#' @return data.table one row with best stats across grid search
#' @export
getBestParameterStats_constraintFDR <- function(complex_feature_grid_filtered,complex_hypotheses,n_subsets=3,FDR_cutoff=0.1){
  # set parameters
  complex_ids <- unique(complex_hypotheses$complex_id)
  targets <- complex_ids[grep("DECOY",complex_ids,invert=TRUE)]
  decoys <- complex_ids[grep("DECOY",complex_ids,invert=FALSE)]
  n_hypotheses <- length(targets)
  subset_size <- floor(n_hypotheses/n_subsets)

  # create random sampling index
  set.seed(123)
  sampled_idx <- sample(seq_len(n_hypotheses), size=n_hypotheses, replace = FALSE)

  # create n subsets of hypotheses
  for(i in seq_len(n_subsets)) {
    start=((i-1)*subset_size)+1
    end=i*subset_size
    idx=sampled_idx[start:end]
    t=targets[idx]
    d=decoys[idx]
    hypotheses=c(t,d)
    name=paste0("subset_",i)
    assign(name,hypotheses)
  }

  # subset results for each n-1 set and estimate FDR
  i_seq=combn(seq_len(n_subsets),n_subsets-1)
  for(i in seq_len(n_subsets)) {
    test_idx=i_seq[,i]
    hypotheses=c()
    for(j in seq_len(length(test_idx))) {
      hypotheses=c(hypotheses,eval(parse(text = paste0("subset_",test_idx[j]))))
    }
    res=lapply(complex_feature_grid_filtered,function(x){subset(x,complex_id %in% hypotheses)})
    stats=estimateGridSearchDecoyFDR(complex_features_list=res)
    stats$P[which(stats$FDR>FDR_cutoff)]=0
    name=paste0("stats_",i)
    assign(name,stats)
  }

  # calculate average P for each parameter setwd
  combi_stats=stats_1
  for(i in seq_len(n_subsets)) {
    if(i>1){
      combi_stats=merge(combi_stats,eval(parse(text = paste0("stats_",i))),
        by=c("corr","window","rt_height","smoothing_length","peak_corr_cutoff","completeness_cutoff","n_subunits_cutoff"),
        suffixes=c("",i))
    }
  }
  p_column_idx=grep("^P",names(combi_stats))
  p_column_names=names(combi_stats)[p_column_idx]
  fdr_column_idx=grep("^FDR",names(combi_stats))
  fdr_column_names=names(combi_stats)[fdr_column_idx]
  combi_stats[,`:=`(meanP = rowMeans(.SD)),by=.I,.SDcols=p_column_names]
  combi_stats[,`:=`(meanFDR = rowMeans(.SD)),by=.I,.SDcols=fdr_column_names]

  combi_stats <- combi_stats[order(-meanP)]

  combi_stats[1]
}
