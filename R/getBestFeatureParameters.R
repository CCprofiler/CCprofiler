#' Perform complex feature grid search filter
#' @description Perform quality filtering of the feature finding result with different
#' cutoffs (grid serch). All possible combinations of the specified parameters are tested.
#' @param grid_search_results List containing result tables from a feature finding grid search.
#' @import data.table
#' @param peak_corr_cutoffs Numeric, vector of within peak correlation_cutoff values to test.
#' Between -1 and 1. Default is c(0.5,0.75,0.9).
#' @param feature_completeness_cutoffs Numeric, vector of feature completeness cutoffs to test.
#' Between 0 and 1. Default is c(0,0.5,1).
#' @param hypothesis_completeness_cutoffs Numeric, vector of hypothesis completeness cutoffs to test.
#' Between 0 and 1. Default is c(0.5,1).
#' @param n_subunits_cutoffs Positive integer vector of minimum number of subunits per hypothesis cutoffs to test. Default is c(2,3,4).
#' @param monomer_distance_cutoffs Positive numeric, factor of allowed distance to monomer weight to test. Default is c(1,2).
#' @param remove_decoys Logical, whether to remove the decoys from the result.
#' Default=\code{FALSE}.
#' @return List of search result tables for every possible parameter combination.
#' The result tables contain additional columns specifying the parameters.
#' @return List with stats
#' @export
#' @examples
#' #------------------------
#' ## Complex level
#' #------------------------
#' 
#' ## Load example data into list to simulate grid search results
#' complexFeaturesGrid <- list(exampleComplexFeatures)
#' 
#' ## Perform the filter grid search
#' complexFeaturesGridFiltered <- filterGridSearchResults(complexFeaturesGrid,
#'                                                        peak_corr_cutoffs = c(0.5,0.75,0.9),
#'                                                        feature_completeness_cutoffs = c(0,0.5,1),
#'                                                        hypothesis_completeness_cutoffs = c(0.5,1),
#'                                                        n_subunits_cutoffs = c(2,3,4),
#'                                                        monomer_distance_cutoffs = c(1,2),
#'                                                        remove_decoys=FALSE)
#' #------------------------
#' ## Protein level
#' #------------------------
#' 
#' ## Load example data into list to simulate grid search results
#' proteinFeaturesGrid <- list(exampleProteinFeatures)
#' 
#' ## Perform the filter grid search
#' proteinFeaturesGridFiltered <- filterGridSearchResults(proteinFeaturesGrid,
#'                                                        peak_corr_cutoffs = c(0.5,0.75,0.9),
#'                                                        feature_completeness_cutoffs = c(0,0.5,1),
#'                                                        hypothesis_completeness_cutoffs = c(0.5,1),
#'                                                        n_subunits_cutoffs =c(2,3,4),
#'                                                        monomer_distance_cutoffs = c(0,1),
#'                                                        remove_decoys=FALSE)
#' 

filterGridSearchResults <- function(grid_search_results,
                                    peak_corr_cutoffs = c(0.5,0.75,0.9),
                                    feature_completeness_cutoffs = c(0,0.5,1),
                                    hypothesis_completeness_cutoffs = c(0.5,1),
                                    n_subunits_cutoffs =c(2,3,4),
                                    monomer_distance_cutoffs = c(1,2),
                                    remove_decoys=FALSE){
  
  .testGridParameter(peak_corr_cutoffs, "peak_corr_cutoffs")
  .testGridParameter(feature_completeness_cutoffs, "feature_completeness_cutoffs")
  .testGridParameter(hypothesis_completeness_cutoffs, "hypothesis_completeness_cutoffs")
  .testGridParameter(n_subunits_cutoffs, "n_subunits_cutoffs")
  .testGridParameter(monomer_distance_cutoffs, "monomer_distance_cutoffs")
  
  parameter_grid <- expand.grid(peak_corr_cutoffs, feature_completeness_cutoffs,
                                hypothesis_completeness_cutoffs, n_subunits_cutoffs,
                                monomer_distance_cutoffs)
  names(parameter_grid) <- c("peak_corr_cutoffs", "feature_completeness_cutoffs",
                             "hypothesis_completeness_cutoffs", "n_subunits_cutoffs",
                             "monomer_distance_cutoffs")
  data <- unlist(lapply(grid_search_results,.helperFilterByParams,
                        params = parameter_grid,
                        remove_decoys = remove_decoys),recursive=FALSE)
  data
}


.helperFilterByParams <- function(data,params,remove_decoys){
  
  res <- apply(params,1,.helperSubset, data = data, remove_decoys = remove_decoys)
  
  res
}


.helperSubset <- function(param,data,remove_decoys){
  
  ## Perform the subset
  x <- filterFeatures(feature_table = data,
                 min_feature_completeness = as.numeric(param["feature_completeness_cutoffs"]),
                 min_hypothesis_completeness = as.numeric(param["hypothesis_completeness_cutoffs"]),
                 min_subunits = as.numeric(param["n_subunits_cutoffs"]),
                 min_peak_corr = as.numeric(param["peak_corr_cutoffs"]),
                 min_monomer_distance_factor = as.numeric(param["monomer_distance_cutoffs"]))
  
  ## Remove decoys if specified
  if(remove_decoys){
    if("protein_id" %in% names(data)){
      res <- subset(res, !grepl("DECOY", protein_id))
    }else if ("complex_id" %in% names(data)){
      res <- subset(res, !grepl("DECOY", complex_id))
    }else{
      stop("Invalid Feature results type. Check the input data!")
    }
  }
  
  ## Add filter criteria columns
  x[,min_feature_completeness := as.numeric(param["feature_completeness_cutoffs"])]
  x[,min_hypothesis_completeness := as.numeric(param["hypothesis_completeness_cutoffs"])]
  x[,min_subunits := as.numeric(param["n_subunits_cutoffs"])]
  x[,min_peak_corr := as.numeric(param["peak_corr_cutoffs"])]
  x[,min_monomer_distance_factor := as.numeric(param["monomer_distance_cutoffs"])]
  
  ## Sort the result
  if("protein_id" %in% names(data)){
    x <- x[order(protein_id,-n_subunits_detected,-area)]
  }else if ("complex_id" %in% names(data)){
    x <- x[order(complex_id,-n_subunits_detected,-area)]
  }else{
    stop("Invalid Feature results type. Check the input data!")
  }
  
  x
}

#' Estimate feature grid search FDR
#' @description Estimate feature grid search FDR by decoy counting.
#' @param complex_features_list data.table containing filtered complex feature results.
#' @param grid_search_params Character vector of column names to report with the statistics for the dataset.
#' Should contain all parameters that are of interest in the grid search.
#' @return A data.table with FDR estimations for every parameter combination specified.
#' @export
#' @examples 
#' ## Generate example data
#' complexFeaturesGrid <- list(exampleComplexFeatures)
#' 
#' complexFeaturesGridFiltered <- filterGridSearchResults(complexFeaturesGrid,
#'                                                        peak_corr_cutoffs = c(0.5,0.75,0.9),
#'                                                        feature_completeness_cutoffs = c(0,0.5,1),
#'                                                        hypothesis_completeness_cutoffs = c(0.5,1),
#'                                                        n_subunits_cutoffs =c(2,3,4),
#'                                                        monomer_distance_cutoffs = c(1,2),
#'                                                        remove_decoys=FALSE)
#' 
#' ## Calculate the FDR statistic for every grid search result
#' # Here we only performed grid search of different filter cutoffs, so only those are
#' # specified as grid_search_params
#' 
#' gridStats <- estimateGridSearchDecoyFDR(complexFeaturesGridFiltered,
#'                                         grid_search_params =c("min_feature_completeness",
#'                                                               "min_hypothesis_completeness",
#'                                                               "min_subunits",
#'                                                               "min_peak_corr",
#'                                                               "min_monomer_distance_factor"))
#' ## Inspect the oputput                                                              
#' gridStats

estimateGridSearchDecoyFDR<- function(complex_features_list, 
                                      grid_search_params =c("corr", "window",
                                                            "rt_height", "smoothing_length",
                                                            "min_feature_completeness",
                                                            "min_hypothesis_completeness",
                                                            "min_subunits",
                                                            "min_peak_corr", 
                                                            "min_monomer_distance_factor")){
  
  x=lapply(complex_features_list,estimateDecoyFDR, 
           grid_search_params = grid_search_params,
           verbose = FALSE)
  
  x_names = names(x[[1]])
  y = as.data.table(t(setDT(x)))[,lapply(.SD,unlist)]
  names(y)=x_names
  y[]
}

#' Get the best performing ParameterStats
#' @description Pick the parameter set in a grid search stats output that performs
#' best in terms of identifications while staying within a specified FDR cutoff.
#' @param grid_search_stats Table of grid search statistics.
#' (obtained from \code{\link{estimateGridSearchDecoyFDR}}).
#' @param FDR Numeric, maximum FDR that should be considered. Default = \code{0.1}.
#' @return data.table one row with best stats across grid search.
#' @export
#' @examples
#' 
#' ## NOT RUN
#' gridStats # see function \code{\link{estimateGridSearchDecoyFDR}} to see how to generate this object.
#' 
#' ## Pick best parameter set
#' getBestParameterStats(gridStats)
#' 

getBestParameterStats <- function(grid_search_stats,
                                  FDR = 0.1){
  env<-environment()
  grid_search_stats <- subset(grid_search_stats,FDR<=get('FDR',env))
  if("P" %in% names(grid_search_stats)) {
    grid_search_stats <- grid_search_stats[order(-P,-TP,FDR)]
  } else {
    grid_search_stats <- grid_search_stats[order(-TP,FDR)]
  }
  grid_search_stats[1]
}


#' Plot FDR gridsearch
#' @description Plot the result of a grid search depending on a specified parameter.
#' @param grid_search_stats Table of grid search statistics 
#' (obtained from \code{\link{estimateGridSearchDecoyFDR}}).
#' @param level Character string, either 'complex' or 'protein'.
#' Specifies which feature finding was performed. Defaults to 'complex'.
#' @param id_level Character string, either 'TP' or 'P'. 
#' Plot with true-positive numbers or all positives as y axis. Defaults to 'TP'
#' @param FDR_cutoff Numeric, the cutoff for the FDR (indicated by a vertical line in the plot).
#' Defaults to \code{0.1}.
#' @param colour_parameter Character string, Which parameter to color. Defaults to 'completeness_cutoff'
#' @param PDF Logical, wether to save the plot as a PDF file in working directory.
#'  Defaults to \code{TRUE}.
#' @param name Character string, filename of the PDF output.
#' @return Either a plot to the R console or a PDF file in the working directory.
#' @export
#' @examples 
#' 
#' ## NOT RUN
#' gridStats # see function \code{\link{estimateGridSearchDecoyFDR}} to see how to generate this object.
#' 
#' ## Plot the result of the grid search depending on the within feature correlation
#' plotIdFDRspace(gridStats, PDF = F, colour_parameter = "min_peak_corr")
#'

plotIdFDRspace <- function(grid_search_stats,
                           level = "complex",
                           id_level = "TP",
                           FDR_cutoff = 0.1,
                           colour_parameter = "min_feature_completeness",
                           PDF = TRUE,
                           name = "ID_FDR_plot.pdf"){
  #@TODO mark best parameter combination and put these parameters in subtitle
  if(level=="complex"){
    sep=100
  } else if (level=="protein"){
    sep=1000
  }
  bestStats <- getBestParameterStats(grid_search_stats,FDR=FDR_cutoff)
  
  sel_best <- which(sapply(1:nrow(grid_search_stats), function(i) identical(grid_search_stats[i,], bestStats)))
  
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


#' Get the best performing parameter data
#' @description Pick the feature-finding result with the parameter set in a grid search output
#' that performs best in terms of identifications while staying within a specified FDR cutoff.
#' @param complex_features_list data.table containing filtered complex feature results.
#' @param FDR Numeric, maximum FDR that should be considered, default = 0.1.
#' @param grid_search_params Character vector of column names to report with the statistics 
#' for the dataset. Should contain all parameters that are of interest in the grid search.
#' @return data.table with features for best parameter set.
#' @export
#' @examples
#' ## Generate example data
#' complexFeaturesGrid <- list(exampleComplexFeatures)
#' 
#' complexFeaturesGridFiltered <- filterGridSearchResults(complexFeaturesGrid,
#'                                                        peak_corr_cutoffs = c(0.5,0.75,0.9),
#'                                                        feature_completeness_cutoffs = c(0,0.5,1),
#'                                                        hypothesis_completeness_cutoffs = c(0.5,1),
#'                                                        n_subunits_cutoffs =c(2,3,4),
#'                                                        monomer_distance_cutoffs = c(1,2),
#'                                                        remove_decoys=FALSE)
#' 
#' ## Extract the result table with the best parameter set
#' bestComplexFeatures <- getBestParameterData(complexFeaturesGridFiltered,
#'                                             grid_search_params = c("min_feature_completeness",
#'                                                                    "min_hypothesis_completeness",
#'                                                                    "min_subunits",
#'                                                                    "min_peak_corr",
#'                                                                    "min_monomer_distance_factor"))
#' head(bestComplexFeatures, n = 2)
#' 

getBestParameterData <- function(complex_features_list,
                                 FDR = 0.1,
                                 grid_search_params =c("corr", "window",
                                                       "rt_height", "smoothing_length",
                                                       "min_feature_completeness",
                                                       "min_hypothesis_completeness",
                                                       "min_subunits",
                                                       "min_peak_corr", 
                                                       "min_monomer_distance_factor")){
  
  grid_stats <- estimateGridSearchDecoyFDR(complex_features_list, 
                                           grid_search_params = grid_search_params)
  
  bestStats <- getBestParameterStats(grid_stats,FDR=FDR)
  
  sel_params <- which(sapply(1:nrow(grid_stats), function(i) identical(grid_stats[i,], bestStats)))
  
  complex_features_list[[sel_params]]
}


#' getBestParameterStats_constraintFDR
#' @description getBestParameterStats_constraintFDR. Do x-fold cross-FDR estimation.
#' Works only for complexFeatures at the moment.
#' @param complex_feature_grid_filtered list with all features from grid search.
#' @param complex_hypotheses data.table with complex hypotheses including decoys.
#' @param n_subsets Numeric, number of subsets used for cross-calculation of FDR, default = 3.
#' @param FDR_cutoff Numeric, maximum FDR that should be considered, default = 0.1.
#' @return data.table one row with best stats across grid search.
#' @export
#' @example 
#' ## Generate example data
#' complexFeaturesGrid <- list(exampleComplexFeatures)
#' complexHypotheses <- exampleComplexHypotheses
#' complexFeaturesGridFiltered <- filterGridSearchResults(complexFeaturesGrid,
#'                                                        peak_corr_cutoffs = c(0.5,0.75,0.9),
#'                                                        feature_completeness_cutoffs = c(0,0.5,1),
#'                                                        hypothesis_completeness_cutoffs = c(0.5,1),
#'                                                        n_subunits_cutoffs =c(2,3,4),
#'                                                        monomer_distance_cutoffs = c(1,2),
#'                                                        remove_decoys=FALSE)
#'  
#' ## Pick the best parameter set with a cross checked FDR constraint of 0.1                                                                                                            
#' getBestParameterStats_constraintFDR(complexFeaturesGridFiltered, 
#'                                     complexHypotheses,
#'                                     FDR = 0.1,
#'                                     grid_search_params =c("min_feature_completeness",
#'                                                           "min_hypothesis_completeness",
#'                                                           "min_subunits",
#'                                                           "min_peak_corr", 
#'                                                           "min_monomer_distance_factor"))
#'                                     
                                    
getBestParameterStats_constraintFDR <- function(complex_feature_grid_filtered,
                                                complex_hypotheses,
                                                n_subsets=3,
                                                FDR_cutoff=0.1,
                                                grid_search_params =c("corr", "window",
                                                                      "rt_height", "smoothing_length",
                                                                      "min_feature_completeness",
                                                                      "min_hypothesis_completeness",
                                                                      "min_subunits",
                                                                      "min_peak_corr", 
                                                                      "min_monomer_distance_factor")){
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
    stats=estimateGridSearchDecoyFDR(complex_features_list=res, 
                                     grid_search_params = grid_search_params)
    stats$P[which(stats$FDR>FDR_cutoff)]=0
    name=paste0("stats_",i)
    assign(name,stats)
  }
  
  # calculate average P for each parameter setwd
  combi_stats=stats_1
  for(i in seq_len(n_subsets)) {
    if(i>1){
      combi_stats=merge(combi_stats,eval(parse(text = paste0("stats_",i))),
                        by=grid_search_params,
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
