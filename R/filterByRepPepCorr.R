#' Filter by sibling peptide correlation
#' @description Filter peptides in traces object based on
#' Replicate Peptide Correlation (RepPepCorr). This function can be called if 
#' one has replicate measurements of one experimental Condition. Trustworthy peptides should
#' show robust elution behaviour over these replicates.
#' @param traces object of class tracesList.
#' @param design_matrix data.table, design matrix describing the architecture of the tracesList object.
#' @param condition_column Character string, The rowname of the design matrix containing
#' the experimental conditions.
#' @param stepwise_filter Logical, wether to iteratively discard replicates until the defined cutoff
#' is met. E.g. if the avg. RepPepCorr is below the cutoff, discard the worst replicate
#' and check again. If the cutoff can't be met with 2 replicates the peptide is discarded.
#' @param repPepCorr_cutoff Numeric, the cutoff below which traces are removed.
#' @return tracesList object containing the integrated traces
#' @export

filterByRepPepCorr <- function(traces,
                               design_matrix,
                               condition_column = NULL,
                               stepwise_filter = TRUE,
                               repPepCorr_cutoff = 0.5){
  
  .tracesListTest(traces)
  
  if(!is.null(condition_column)){
    compare_col <- which(names(design_matrix) == condition_column)
    stopifnot(length(compare_col) == 1)
    conditions <- unique(design_matrix[[condition_column]])
  }else{
    design_matrix$Condition <- "all Conditions"
    conditions <- "all Conditions"
    condition_column <- "Condition"
  }
  if (!("RepPepCorr" %in% names(traces$trace_annotation))){
    message("Replicate peptide correlation not yet calculated for this dataset\nCalculating RepPepCorr...")
    traces <- calculateRepPepCorr(traces, design_matrix, compare_within = condition_column, add = T)
  } else{
    message("Replicate peptide correlation values found...")
  }
  
  message(paste("Applying RepPepCorr cutoff of ", repPepCorr_cutoff))
  if(stepwise_filter){
    #Filter each sample separately
    traces_res <- lapply(traces, function(tr){
      ids <- tr$trace_annotation[RepPepCorr >= repPepCorr_cutoff, id]
      subset(tr, ids)
    })
  }else{
    traces_res <- list()
    #Filter with the mean RepPeporr in every condition
    for(condition in conditions){
      samplesInCond <- design_matrix[get(compare_within) == condition, Sample_name]
      corrs <- do.call(rbind, lapply(traces[samplesInCond], "[[", "trace_annotation"))
      meancorrs <- corrs[,.(meanRPC = mean(RepPepCorr)), by = id]
      ids <- meancorrs[meanRPC >= repPepCorr_cutoff, id]
      for(sample in samplesInCond){
        traces_res[[sample]] <- subset(traces[[sample]], ids)
      }
    }
  }
  class(traces_res) <- "tracesList"
  .tracesListTest(traces_res)
  return(traces_res)
  
}


