#' Estimate statistics based on decoy model.
#' @description Estimate statistics based on decoy model.
#' @param complex_features data.table containing filtered complex feature results.
#' @return List with stats
#' @export
estimateDecoyFDR <- function(complex_features,grid_search_list=FALSE){
  complexes <- complex_features$complex_id
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
#' @param complex_features data.table containing filtered complex feature results.
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
  name=paste(params,collapse="_")
  name=gsub("\\.","",name)
  name=gsub(" ","",name)
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
#' @return List with stats
#' @export
plotIdFDRspace <- function(grid_search_stats,FDR_cutoff=0.1,colour_parameter="completeness_cutoff",PDF=TRUE,name="ID_FDR_plot.pdf"){
  #@TODO mark best parameter combination and put these parameters in subtitle
  if(PDF){pdf(name)}
  pl <- ggplot(data=grid_search_stats,aes(y=TP,x=FDR,colour=get(colour_parameter))) +
    geom_point() +
    scale_x_continuous(breaks=seq(0,1,0.1),limits=c(0,1),minor_breaks=NULL) +
    scale_y_continuous(breaks=seq(0,ceiling(max(grid_search_stats$TP)/100)*100,100),limits=c(0,ceiling(max(grid_search_stats$TP)/100)*100),minor_breaks=NULL) +
    scale_colour_continuous(guide = guide_legend(title = eval(colour_parameter))) +
    geom_vline(xintercept=FDR_cutoff,colour="red",linetype=2) +
    theme_bw()
  print(pl)
  if(PDF){dev.off()}
}
