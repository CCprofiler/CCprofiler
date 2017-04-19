#' Estimate statistics based on decoy model.
#' @description Estimate statistics based on decoy model.
#' @param complex_features data.table containing filtered complex feature results.
#' @return List with stats
#' @export
estimateDecoyFDR <- function(complex_features,grid_search_list=FALSE){
  if("complex_id" %in% names(complex_features)){
    complexes <- complex_features$complex_id
  } else if ("protein_id" %in% names(complex_features)) {
    complexes <- complex_features$protein_id
  } else {
    message("not a valid search result")
  }
  decoys = complexes[grep("DECOY",complexes)]
  targets = complexes[!complexes %in% decoys]
  if(!length(decoys)>0){
    message("No decoy features found. Please check if decoys were available in the complex hypotheses.")
    FDR = 0
  } else {
    # estimate FDR based on detected decoys
    FDR = length(decoys)/length(targets)
  }
  # estimate number of true positive complexes by correcting with the estimated FDR
  P = length(targets)
  TP = length(targets)*(1-FDR)
  if(grid_search_list){
    list(FDR=FDR,TP=TP,P=P,corr=unique(complex_features$corr),window=unique(complex_features$window),rt_height=unique(complex_features$rt_height),smoothing_length=unique(complex_features$smoothing_length),peak_corr_cutoff=unique(complex_features$peak_corr_cutoff),completeness_cutoff=unique(complex_features$completeness_cutoff),n_subunits_cutoff=unique(complex_features$n_subunits_cutoff))
  } else {
    list(FDR=FDR,TP=TP,P=P)
  }
}

#' Perform complex feature grid search
#' @description Perform complex feature grid search.
#' @param traces traces object of type protein
#' @param calibration list of two functions for calibration
#' @param corrs numeric vector
#' @param windows numeric vector
#' @param smoothing numeric vector
#' @param rt_heights numeric vector
#' @param parallelized logical default=TRUE
#' @param n_cores numeric number of cores to use if parallelized (default=1)
#' @return List with stats
#' @export
performComplexGridSearch <- function(traces,
                                    complex_hypothesis,
                                    calibration,
                                    corrs = c(0.5,0.75,0.9,0.95),
                                    windows = c(8,10,12),
                                    smoothing = c(7,9,11),
                                    rt_heights = c(3,4,5),
                                    parallelized = TRUE,
                                    n_cores=1
                                    ){
  parameter_grid <- expand.grid(corrs,windows,smoothing,rt_heights)
  names(parameter_grid) <- c("corr","window","smoothing","rt_height")
  cl <- snow::makeCluster(n_cores)
  # setting a seed is absolutely crutial to ensure reproducible results!!!!!!!!!!!!!!!!!!!
  clusterSetRNGStream(cl,123)
  doSNOW::registerDoSNOW(cl)
  clusterEvalQ(cl,library(SECprofiler,data.table))
  data <- list()
  data <- c(data,parRapply(cl,parameter_grid,FUN=runGridComplexFeatureFinding,protTraces=traces,calibration=calibration,complex_hypothesis=complex_hypothesis))
  stopCluster(cl)
  data
}

runGridComplexFeatureFinding <- function(params,protTraces,calibration,complex_hypothesis) {
  res = findComplexFeatures(traces=protTraces,
                            complex_hypothesis = complex_hypothesis,
                            calibration = calibration,
                            corr_cutoff = as.numeric(params["corr"]),
                            window_size = as.numeric(params["window"]),
                            parallelized = FALSE,
                            collapse_method="apex_network",
                            perturb_cutoff="1%",
                            rt_height=as.numeric(params["rt_height"]),
                            smoothing_length=as.numeric(params["smoothing"])
                            )
  res[,corr:=as.numeric(params["corr"])]
  res[,window:=as.numeric(params["window"])]
  res[,rt_height:=as.numeric(params["rt_height"])]
  res[,smoothing_length:=as.numeric(params["smoothing"])]
  #saveRDS(res,file=paste0("complexFeatures/complexFeatures_",name,".rda"))
  res[]
}

#' Perform complex feature grid search filter
#' @description Perform complex feature grid search filter.
#' @param complex_features data.table containing filtered complex feature results.
#' @return List with stats
#' @export
filterGridSearchResults <- function(grid_search_results,
                                    peak_corr_cutoffs = c(0.5,0.75,0.9),
                                    completeness_cutoffs = c(0,0.5,1),
                                    n_subunits_cutoffs =c(2,3,4),
                                    remove_decoys=FALSE
                                    ) {
  parameter_grid <- expand.grid(peak_corr_cutoffs,completeness_cutoffs,n_subunits_cutoffs)
  names(parameter_grid) <- c("peak_corr_cutoffs","completeness_cutoffs","n_subunits_cutoffs")
  data <- list()
  data <- c(data,unlist(lapply(grid_search_results,helperFilterByParams,params=parameter_grid,remove_decoys=remove_decoys),recursive=FALSE))
  data
}

helperFilterByParams <- function(data,params,remove_decoys){
  res <- list()
  res <- c(res,apply(params,1,helperSubset,data=data,remove_decoys=remove_decoys))
  res
}

helperSubset <- function(param,data,remove_decoys){
  x <- subset(data,((peak_corr >= as.numeric(param["peak_corr_cutoffs"])) & (completeness >= as.numeric(param["completeness_cutoffs"])) & (n_subunits_detected >= as.numeric(param["n_subunits_cutoffs"]))))
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
  x[]
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

#' Plot FDR gridsearch
#' @description Perform complex feature grid search.
#' @param complex_features data.table containing filtered complex feature results.
#' @param level character sting either complex or protein
#' @return List with stats
#' @export
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
