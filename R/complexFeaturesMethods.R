
resultsToTable <- function(swf){
  res <- swf$sw.results
  res.list <- lapply(seq(1:length(res)), function(i){
    features <- res[[i]]$features
    features
  })
  res <- do.call("rbind", res.list)        
  res
}

subsetComplexFeatures <- function(res,complex_ids=NULL,min_completeness=NULL){
  if(!is.null(complex_ids)){
    res <- subset(res,complex_id %in% as.character(complex_ids))
  }
  if(!is.null(min_completeness)){
    res <- subset(res,completeness >= min_completeness)
  }
  res
}

getBestComplexFeature <- function(swf){
  res <- swf$sw.results
  res.list <- lapply(seq(1:length(res)), function(i){
    features <- res[[i]]$features
    features[1] 
  })
  res <- do.call("rbind", res.list)
  res
} 


