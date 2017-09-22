calculatePairwiseRepPepCorr <- function(traces, comparisons = NULL){
  
  # Check input type
  if (any(sapply(traces, "[[", "trace_type") != "peptide")){
    stop("Replicate peptide correlation can only be calculated on traces of type peptide")
  }
  for(i in 1:length(traces)){
    setkey(traces[[i]]$traces, "id")
  }
  if(is.null(comparisons)){
    corrPairs <- combn(names(traces), m = 2)
  }else{
    corrPairs <- comparisons
  }
  
  corrPairNames <- apply(corrPairs, 2, paste, collapse = "-")
  intMats <- lapply(traces, getIntensityMatrix)
  allIds <- unique(do.call("c", lapply(intMats, rownames)))
  res <- data.table(id = allIds)
  
  cross_corr <- apply(corrPairs, 2, function(pair){
    # print(pair)
    tr1 <- intMats[[pair[1]]]
    tr2 <- intMats[[pair[2]]]
    ids1 <- rownames(tr1)
    ids2 <- rownames(tr2)
    common_ids <- intersect(ids1, ids2)
    
    t1 <- tr1[common_ids,]
    t2 <-  tr2[common_ids,]
    # sapply(1:nrow(t1), function(i) cor(t1[i,], t2[i,]))
    cA <- t1 - rowMeans(t1)
    cB <- t2 - rowMeans(t2)
    sA <- sqrt(rowMeans(cA^2))
    sB <- sqrt(rowMeans(cB^2))
    
    rowMeans(cA * cB) / (sA * sB)
  })
  cross_corr <- lapply(cross_corr, function(x) data.table(id = names(x), corr = x))
  for(i in 1: length(cross_corr)){
    res <- merge(res, cross_corr[[i]], by = "id", all.x = TRUE, suffixes = c(i-1, i))
  }
  names(res) <- c("id", corrPairNames)
  resLong <- melt(res, id.vars = "id", variable.name = "Sample_pair", value.name = "RepPepCorr")
  ann <- as.data.table(cbind(corrPairNames, t(corrPairs)))
  names(ann) <- c("Sample_pair", "Sample1", "Sample2")
  resLong <- merge(resLong, ann, by = "Sample_pair", all.x = TRUE)
  setkey(resLong, "id")
  return(resLong)  
}


integrateReplicates <- function(traces,
                                design_matrix,
                                integrate_within = NULL,
                                filter_by_RepPepCorr = TRUE,
                                repPepCorr_cutoff = 0.5){
  
  
  
  if(!is.null(integrate_within)){
    compare_col <- which(names(design_matrix) == integrate_within)
    stopifnot(length(compare_col) == 1)
    conditions <- unique(design_matrix[[integrate_within]])
  }else{
    design_matrix$Condition <- "all Conditions"
    conditions <- "all Conditions"
    integrate_within <- "Condition"
  }
  
  traces_res <- list()
  pairCorr <- calculatePairwiseRepPepCorr(traces = traces)
  for(condition in conditions){
    message(paste("Integrating rplicates within:", condition))
    sample_names <- design_matrix[get(integrate_within) == condition, Sample_name]
    
    if(filter_by_RepPepCorr){
      
      pc <- pairCorr[Sample1 %in% sample_names & Sample2 %in% sample_names & !is.na(RepPepCorr)]
      pc[, meanRepPepCorr := mean(RepPepCorr), by = id]
      pcPass <- pc[meanRepPepCorr >= repPepCorr_cutoff]
      
      removeWorst <- function(tabl){
        a <- melt(tabl, measure.vars = c("Sample1", "Sample2"))[, variable :=NULL]
        setkey(a, "id")
        a[, repPepCorrSum := sum(RepPepCorr), by = .(id, value)]
        a <- a[a[,repPepCorrSum != min(repPepCorrSum), by = id]$V1]
        a <- unique(a[,.(id, value)])
        
        res <- apply(a, 1, function(x) tabl[id == x[1] & Sample1 != x[2] & Sample2 != x[2]]) 
        res <- do.call("rbind", res)
        res[, meanRepPepCorr := mean(RepPepCorr), by = id]
        res
      }
      i <- length(sample_names)
      while(i >2 & nrow(pc) > 0){
        # Remove the worst replicate
        pc <- pc[meanRepPepCorr < repPepCorr_cutoff]
        pc <- removeWorst(pc)
        pcPass <- rbind(pcPass, pc[meanRepPepCorr >= repPepCorr_cutoff])
        i <- i-1
      }
      setkey(pcPass, "id")
      
      for(sample_name in sample_names){
        subset <- unique(pcPass[Sample1 == sample_name | Sample2 == sample_name, id])
        traces[[sample_name]] <- subset(traces[[sample_name]], trace_subset_ids = )
      }
      
    }
    
    traces_res[[condition]] <- sumTraceIntensities(traces = traces[sample_names], aggr_fun = "mean")
  }
  class(traces_res) <- "tracesList"
  .tracesListTest(traces_res)
  return(traces_res)
}

sumTraceIntensities <- function(traces, min_present = 1, aggr_fun = c("sum", "mean")){
  
  for(i in 1:length(traces)){
    .tracesTest(traces[[i]], type = traces[[1]]$trace_type)
  }
  aggr_fun <- match.arg(aggr_fun)
  
  ## Create a merged annotation
  annComb <- do.call("rbind", lapply(traces, function(x) x$trace_annotation))
  setkey(annComb, "id")
  if(any(names(annComb) == "SibPepCorr")){
    annComb[, meanSibPepCorr := mean(SibPepCorr), by= id][, SibPepCorr := NULL]
  }
  if(any(names(annComb) == "RepPepCorr")){
    annComb[, meanRepPepCorr := mean(RepPepCorr), by= id][, RepPepCorr := NULL]
  }
  annComb[, detectedIn := .N, by = id]
  annComb <- unique(annComb, by = "id")
  
  ## Calculate Mean of intensities
  traceInt <- do.call("rbind", lapply(traces, "[[", "traces"))
  setkey(traceInt, "id")
  if(aggr_fun == "mean"){
    traceIntSum <- traceInt[, lapply(.SD,function(x) sum(x)/.N), by = id]
  }else if(aggr_fun == "sum"){
    traceIntSum <- traceInt[, lapply(.SD, sum), by = id]
  }
  setcolorder(traceIntSum, c(names(traceIntSum)[-1], "id"))
  
  ## Create a merged fraction annotation
  fracAnnComb <- rbindlist(lapply(traces, "[[", "fraction_annotation"), fill = TRUE)
  setkey(fracAnnComb, "id")
  if(ncol(fracAnnComb) == 1){
    fracAnnComb <- unique(fracAnnComb)
  }else{
    fracAnnComb <- fracAnnComb[, lapply(.SD, function(x){
      if(all(x == x[1],na.rm = T)){
        x[1]
      } else{
        paste0(x, collapse = ",")
      }
    }), by = id]
  }
   
  ## Assemble the new traces object
  resTraces <- list(traces = traceIntSum,
                    trace_type = traces[[1]]$trace_type,
                    trace_annotation = annComb,
                    fraction_annotation = fracAnnComb)
  class(resTraces) <- "traces"
  .tracesTest(resTraces)
  
  return(resTraces)
}
