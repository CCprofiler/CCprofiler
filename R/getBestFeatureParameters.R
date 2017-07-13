#' Perform complex feature grid search filter
#' @description Perform complex feature grid search filter.
#' @param grid_search_results data.table containing filtered complex feature results.
#' @return List with stats
#' @export
#' 

filterGridSearchResults <- function(grid_search_results,
                                    peak_corr_cutoffs = c(0.5,0.75,0.9),
                                    completeness_cutoffs = c(0,0.5,1),
                                    n_subunits_cutoffs =c(2,3,4),
                                    remove_decoys=FALSE
) {
  parameter_grid <- expand.grid(peak_corr_cutoffs,completeness_cutoffs,n_subunits_cutoffs)
  names(parameter_grid) <- c("peak_corr_cutoffs","completeness_cutoffs","n_subunits_cutoffs")
  data <- list()
  data <- c(data,unlist(lapply(grid_search_results,helperFilterByParams,
                               params = parameter_grid,
                               remove_decoys = remove_decoys),recursive=FALSE))
  data
}

helperFilterByParams <- function(data,params,remove_decoys){
  res <- list()
  res <- c(res,apply(params,1,helperSubset,data=data,remove_decoys=remove_decoys))
  res
}

helperSubset <- function(param,data,remove_decoys){
  x <- subset(data,((peak_corr >= as.numeric(param["peak_corr_cutoffs"])) & (n_subunits_detected >= as.numeric(param["n_subunits_cutoffs"]))))
  x <- x[order(complex_id,-n_subunits_detected,-sw_score,-area,mw_diff)]
  if("complex_id" %in% names(x)){
    allowed_ids <- x[completeness>=as.numeric(param["completeness_cutoffs"]), unique(complex_id)]
    x <- x[complex_id %in% allowed_ids]
  } else if ("protein_id" %in% names(x)) {
    allowed_ids <- x[completeness>=as.numeric(param["completeness_cutoffs"]), unique(protein_id)]
    x <- x[protein_id %in% allowed_ids]
  } else {
    message("not a valid search result")
  }
  if (remove_decoys) {
    if("complex_id" %in% names(x)){
      x <- x[!(grep("DECOY",complex_id))]
    } else if ("protein_id" %in% names(x)) {
      x <- x[!(grep("DECOY",protein_id))]
    } else {
      message("not a valid search result")
    }
  }
  x[,peak_corr_cutoff:=as.numeric(param["peak_corr_cutoffs"])]
  x[,completeness_cutoff:=as.numeric(param["completeness_cutoffs"])]
  x[,n_subunits_cutoff:=as.numeric(param["n_subunits_cutoffs"])]
  x <- x[order(complex_id,-n_subunits_detected,-sw_score,-area,mw_diff)]
  x[]
}


#' Plot FDR gridsearch
#' @description Perform complex feature grid search.
#' @param complex_features data.table containing filtered complex feature results.
#' @param level character sting either complex or protein
#' @return List with stats
#' @export
#' 

plotIdFDRspace <- function(grid_search_stats,level="complex",id_level="TP",FDR_cutoff=0.1,colour_parameter="completeness_cutoff",PDF=TRUE,name="ID_FDR_plot.pdf"){
  #@TODO mark best parameter combination and put these parameters in subtitle
  if(level=="complex"){
    sep=100
  } else if (level=="protein"){
    sep=1000
  }
  BestStats <- getBestParameterStats(grid_search_stats,FDR=FDR_cutoff)
  sel_best <- which((grid_search_stats$peak_corr_cutoff==BestStats$peak_corr_cutoff) &
                      (grid_search_stats$corr==BestStats$corr) &
                      (grid_search_stats$window==BestStats$window) &
                      (grid_search_stats$rt_height==BestStats$rt_height) &
                      (grid_search_stats$smoothing_length==BestStats$smoothing_length) &
                      (grid_search_stats$peak_corr_cutoff==BestStats$peak_corr_cutoff) &
                      (grid_search_stats$completeness_cutoff==BestStats$completeness_cutoff) &
                      (grid_search_stats$n_subunits_cutoff==BestStats$n_subunits_cutoff))
  grid_search_stats[,best:=1]
  grid_search_stats$best[sel_best] = 1.5
  if(id_level=="TP"){
    if(PDF){pdf(name)}
    pl <- ggplot(data=grid_search_stats,aes(y=TP,x=FDR,colour=factor(get(colour_parameter)),size=best)) +
      geom_point() +
      scale_x_continuous(breaks=seq(0,1,0.1),limits=c(0,1),minor_breaks=NULL) +
      scale_y_continuous(breaks=seq(0,ceiling(max(grid_search_stats$TP)/100)*100,sep),limits=c(0,ceiling(max(grid_search_stats$TP)/100)*100),minor_breaks=NULL) +
      scale_size(guide = FALSE) +
      scale_colour_hue(guide = guide_legend(title = paste0(eval(colour_parameter),"\n"))) +
      geom_vline(xintercept=FDR_cutoff,colour="red",linetype=2) +
      theme_bw()
    print(pl)
    if(PDF){dev.off()}
  } else if (id_level=="P"){
    if(PDF){pdf(name)}
    pl <- ggplot(data=grid_search_stats,aes(y=P,x=FDR,colour=factor(get(colour_parameter)),size=best)) +
      geom_point() +
      scale_x_continuous(breaks=seq(0,1,0.1),limits=c(0,1),minor_breaks=NULL) +
      scale_y_continuous(breaks=seq(0,ceiling(max(grid_search_stats$P)/100)*100,sep),limits=c(0,ceiling(max(grid_search_stats$P)/100)*100),minor_breaks=NULL) +
      scale_size(guide = FALSE) +
      scale_colour_hue(guide = guide_legend(title = paste0(eval(colour_parameter),"\n"))) +
      geom_vline(xintercept=FDR_cutoff,colour="red",linetype=2) +
      theme_bw()
    print(pl)
    if(PDF){dev.off()}
  } else {
    message("pick TP or P as id_level")
  }
}

#' getBestParameterStats
#' @description getBestParameterStats.
#' @param stats data.table with all stats from grid search
#' @param FDR numeric maximum FDR tat should be considered, default = 0.1
#' @return data.table one row with best stats across grid search
#' @export
getBestParameterStats <- function(stats,FDR=0.1){
  env<-environment()
  stats <- subset(stats,FDR<=get('FDR',env))
  if("P" %in% names(stats)) {
    stats <- stats[order(-P,-TP,FDR)]
  } else {
    stats <- stats[order(-TP,FDR)]
  }
  stats[1]
}

#' getBestParameterData
#' @description getBestParameterData.
#' @param grid_data list with all filtered grid search results
#' @param FDR numeric maximum FDR tat should be considered, default = 0.1
#' @return data.table with complex features for best parameter set
#' @export
getBestParameterData <- function(grid_data,FDR=0.1){
  grid_stats <- estimateGridSearchDecoyFDR(grid_data)
  BestStats <- getBestParameterStats(grid_stats,FDR=FDR)
  sel_params <- which((grid_stats$peak_corr_cutoff==BestStats$peak_corr_cutoff) &
                        (grid_stats$corr==BestStats$corr) &
                        (grid_stats$window==BestStats$window) &
                        (grid_stats$rt_height==BestStats$rt_height) &
                        (grid_stats$smoothing_length==BestStats$smoothing_length) &
                        (grid_stats$peak_corr_cutoff==BestStats$peak_corr_cutoff) &
                        (grid_stats$completeness_cutoff==BestStats$completeness_cutoff) &
                        (grid_stats$n_subunits_cutoff==BestStats$n_subunits_cutoff))
  grid_data[[sel_params]]
}

#' Perform complex feature grid search FDR
#' @description Perform complex feature grid search.
#' @param complex_features data.table containing filtered complex feature results.
#' @return List with stats
#' @export
estimateGridSearchDecoyFDR<- function(complex_features_list){
  x=lapply(complex_features_list,estimateDecoyFDR,grid_search_list=TRUE)
  x_names = names(x[[1]])
  y = as.data.table(t(setDT(x)))[,lapply(.SD,unlist)]
  names(y)=x_names
  y[]
}

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
