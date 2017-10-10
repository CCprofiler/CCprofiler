#' Calculate coelution score for all detected features.
#' @param features data.table with complex or protein features
#' @return data.table with complex or protein features including an extra coelution_score column
#' @export
calculateCoelutionScore <- function(features){
  features[peak_corr < 0]$peak_corr = 0
  #features <- subset(features,peak_corr > 0)
  features$coelution_score <- apply(features,1,function(x){
  peak_corr=as.numeric(x["peak_corr"])
  n_subunits_detected=as.numeric(x["n_subunits_detected"])
  n_subunits_annotated=as.numeric(x["n_subunits_annotated"])
  1-min(1,sum(((1-peak_corr)^((n_subunits_detected:n_subunits_annotated)-1))*((peak_corr)^(n_subunits_annotated-(n_subunits_detected:n_subunits_annotated)))*choose(n_subunits_annotated-1,(n_subunits_detected:n_subunits_annotated)-1)))
  })
  return(features[])
}

#' Calculate q-values from coelution score for all detected features.
#' @param features data.table with complex or protein features
#' @param lambda numeric value between 0 and 1, default = 0.5
#' @param plot logical, default = TRUE
#' @param PDF logical, defalt = FALSE
#' @param name character strimg specifying pdf name, default = "q_value_stats"
#' @importFrom qvalue empPvals
#' @importFrom qvalue qvalue
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
#' @export
calculateQvalue <- function(features,lambda=0.5,plot=TRUE,PDF=FALSE,name="q_value_stats"){
  features[,decoy:=0]
  if ("complex_id" %in% names(features)) {
    features$decoy[grep("DECOY",features$complex_id)] = 1
  } else if ("protein_id" %in% names(features)) {
    features$decoy[grep("DECOY",features$protein_id)] = 1
  } else {
    stop("Not a complex of protein feature table.")
  }
  targets=features[decoy==0,coelution_score]
  decoys=features[decoy==1,coelution_score]
  pvalues <- empPvals(stat=targets,stat0=decoys,pool=TRUE)
  qobj <- qvalue(pvalues,lambda=lambda)
  if(plot){
    if(PDF){
      pdf(paste0(name,".pdf"))
    }
      hist(pvalues,nclass=20)
      plot(qobj)
      print(hist(qobj))
    if(PDF){
      dev.off()
    }
  }
  target_data <- subset(features,decoy==0)
  target_data[,pvalue:=pvalues]
  target_data[,qvalue:=qobj$qvalues]
  return(target_data[])
}

#' qvaluePositivesPlot
#' @param features data.table with complex or protein features including qvalues
#' @param plot logical, default = TRUE
#' @param PDF logical, defalt = FALSE
#' @param name character strimg specifying pdf name, default = "qvaluePositivesPlot"
#' @return data.table with stats
#' @export
qvaluePositivesPlot <- function(features,plot=TRUE,PDF=FALSE,name="qvaluePositivesPlot"){
  maxQvalue <- max(features$qvalue)
  stats <- data.table(qvalue_cutoff=seq(0,maxQvalue,0.005))
  stats[,positives := unlist(lapply(qvalue_cutoff,function(x)nrow(subset(features,qvalue<=x))))]
  stats[,true_positives := unlist(lapply(qvalue_cutoff,function(x){nrow(subset(features,qvalue<=x))*(1-x)}))]
  if(plot){
    if(PDF){
      pdf(paste0(name,".pdf"))
    }
      p <- ggplot(data=stats,aes(x=qvalue_cutoff,y=positives)) +
        geom_point()
      print(p)
      tp <- ggplot(data=stats,aes(x=qvalue_cutoff,y=true_positives)) +
        geom_point()
      print(tp)
    if(PDF){
      dev.off()
    }
  }
  return(stats[])
}

#' plotScoreDistribution
#' @param features data.table with complex or protein features including qvalues
#' @param PDF logical, defalt = FALSE
#' @param name character strimg specifying pdf name, default = "scoreDistribution"
#' @export
plotScoreDistribution <- function(features,PDF=TRUE,name="scoreDistribution"){
  pdf(paste0(name,".pdf"))
   features[,decoy:=0]
   if ("complex_id" %in% names(features)) {
     features$decoy[grep("DECOY",features$complex_id)] = 1
   } else if ("protein_id" %in% names(features)) {
     features$decoy[grep("DECOY",features$protein_id)] = 1
   } else {
     stop("Not a complex of protein feature table.")
   }
    pl <- ggplot(data=features,aes(x=coelution_score,fill=factor(decoy))) + geom_histogram(position="dodge",binwidth=0.05)+
    scale_x_continuous(breaks=seq(0,1,0.1),limits=c(0,1),minor_breaks=NULL)
    print(pl)
  dev.off()
}

#' qvaluePositivesPlotGrid
#' @param featuresGrid list of complex or protein feature data.tables including qvalues
#' @param plot logical, default = TRUE
#' @param PDF logical, defalt = FALSE
#' @param name character strimg specifying pdf name, default = "qvaluePositivesPlotGrid"
#' @return data.table with stats
#' @export
qvaluePositivesPlotGrid <- function(featuresGrid,plot=TRUE,PDF=FALSE,name="qvaluePositivesPlotGrid"){
  stats_list <- lapply(featuresGrid,function(x){
    maxQvalue <- max(x$qvalue)
    stats <- data.table(qvalue_cutoff=seq(0,maxQvalue,0.005))
    stats[,positives := unlist(lapply(qvalue_cutoff,function(y)nrow(subset(x,qvalue<=y))))]
    stats[,true_positives := unlist(lapply(qvalue_cutoff,function(y){nrow(subset(x,qvalue<=y))*(1-y)}))]
    stats[,corr := x$corr[1]]
    stats[,window := x$window[1]]
    stats[,rt_height := x$rt_height[1]]
    stats[,smoothing_length := x$smoothing_length[1]]
    stats[]
  })
  stats_all <- do.call("rbind", stats_list)
  if(plot){
    if(PDF){
      pdf(paste0(name,".pdf"))
    }
      p <- ggplot(data=stats_all,aes(x=qvalue_cutoff,y=positives,colour=factor(corr))) +
        geom_point()
      print(p)
      tp <- ggplot(data=stats_all,aes(x=qvalue_cutoff,y=true_positives,colour=factor(corr))) +
        geom_point()
      print(tp)
    if(PDF){
      dev.off()
    }
  }
  return(stats_all)
}

#' getBestQvalueParameters
#' @param stats data.table output from qvaluePositivesPlotGrid function
#' @param FDR_cutoff numeric between 0 and 1, default = 0.05
#' @return data.table with best stats
#' @export
getBestQvalueParameters <- function(stats,FDR_cutoff=0.05){
  stats_sub <- subset(stats,qvalue_cutoff<=FDR_cutoff)
  stats_sub <- stats_sub[order(-positives)]
  best_stats <- stats_sub[1]
  return(best_stats)
}
