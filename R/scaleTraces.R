
scaleTraces <- function(traces, design_matrix, common_scale_within = "Condition"){
  .tracesListTest(traces)
  max_values <- lapply(traces, function(tr){
    max(getIntensityMatrix(tr))
  })
  conditions <- unique(design_matrix[[common_scale_within]])
  rescale_vals <- lapply(conditions, function(cond){
    samples <- design_matrix[get(common_scale_within) == cond, Sample_name]
    mean(unlist(max_values[samples]))
  })
  names(rescale_vals) <- conditions
  
  res <- lapply(names(traces), function(sample){
    cond <- design_matrix[Sample_name == sample, get(common_scale_within)]
    norm_val <- rescale_vals[[cond]] / max_values[[sample]]
    normalizeTraces(traces[[sample]], 
                    normVals = rep(norm_val, nrow(traces[[sample]]$fraction_annotation)),
                    method = "*", transform = "none")
  })
  names(res) <- names(traces)
  class(res) <- "tracesList"
  .tracesListTest(res)
  return(res)
}
