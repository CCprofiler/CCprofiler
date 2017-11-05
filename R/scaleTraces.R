
#' Normalize tracesList intensities to a common scale
#' @description 
#' @param traces Object of class tracesList.
#' @param design_matrix data.table, A valid design matrix describing the tracesList object.
#' @param common_scale_within Character string, the name of a column in the design matrix.
#' Values within groups in that column are brought to a common scale. Defaults to 'Condition'.
#' @param verbose Logical, wether to print messages to the console.
#' @return Object of class traces with normalized intensities.
#' @export

scaleTraces <- function(traces, design_matrix,
                        common_scale_within = "Condition",
                        scaling_method = c("max", "median", "mean"),
                        verbose = TRUE){
  .tracesListTest(traces)
  scaling_method <- match.arg(scaling_method)
  scale_fun <- match.fun(scaling_method)
  summary_values <- lapply(traces, function(tr){
    m <- getIntensityMatrix(tr)
    scale_fun(m[m > 0])
  })
  if(!is.null(common_scale_within)){
    conditions <- unique(design_matrix[[common_scale_within]])
  } else{
    conditions <- 1
  }
  
  rescale_vals <- lapply(conditions, function(cond){
    if(!is.null(common_scale_within)){
      samples <- design_matrix[get(common_scale_within) == cond, Sample_name]
    }else{
      samples <- design_matrix[,Sample_name]
    }
    mean(unlist(summary_values[samples]))
  })
  names(rescale_vals) <- conditions
  
  res <- lapply(names(traces), function(sample){
    if(!is.null(common_scale_within)){
    cond <- design_matrix[Sample_name == sample, get(common_scale_within)]
    norm_val <- rescale_vals[[cond]] / summary_values[[sample]]
    } else{
      norm_val <- rescale_vals[[1]] / summary_values[[sample]]
    }
    
    if(verbose) message(paste0("Scaling ", sample, " with factor ", norm_val))
    normalizeTraces(traces[[sample]], 
                    normVals = rep(norm_val, nrow(traces[[sample]]$fraction_annotation)),
                    method = "*", transform = "none")
  })
  names(res) <- names(traces)
  class(res) <- "tracesList"
  .tracesListTest(res)
  return(res)
}
