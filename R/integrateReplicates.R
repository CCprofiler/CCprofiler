## Legacy function

#' Generate a consensus traces object from multiple replicate traces
#' @param traces object of class tracesList.
#' @param design_matrix data.table, design matrix describing the architecture of the tracesList object.
#' @param integrate_within Character string, Condition to collapse 
#' (must be a valid column of the design matrix). If \code{NULL}, everything is collapsed.
#' @param filter_by_RepPepCorr Logical, wether to apply a RepPepCorr filter.
#' @param repPepCorr_cutoff Numeric, the cutoff below which traces are removed.
#' @return tracesList object containing the integrated traces
#' @export
integrateReplicates <- function(traces,
                                design_matrix,
                                integrate_within = NULL,
                                filter_by_RepPepCorr = TRUE,
                                repPepCorr_cutoff = 0.5){
  
  .tracesListTest(traces)
  
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
    message(paste("Integrating replicates within:", condition))
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

