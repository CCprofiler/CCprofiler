
#' Normalize tracesList intensities to a common scale
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


#' Plot a comparison of the intensities in different samples
#' @param traces Object of class tracesList.
#' @param PDF Logical, whether to plot to PDF. PDF file is saved in working directory.
#'  Default is \code{FALSE}.
#' @param name Character string with name of the plot, only used if \code{PDF=TRUE}.
#' PDF file is saved under name.pdf. Default is "Traces".
#' @param plot Logical, wether to print or return the plot object
#' @return Violin plot of trace intensities per sample
#' @export

plotGlobalIntensities <- function(traces, plot = T, PDF = F, name = "IntensitySummary.pdf"){
  .tracesListTest(traces)
  allInts <- do.call(rbind, lapply(names(traces), function(x){
    m <- c(getIntensityMatrix(traces[[x]]))
    data.table( Intensity = m[m>0], Sample = x)
  }))
  means <- allInts[, .(Intensity = mean(Intensity)), by = Sample]
  p <- ggplot(allInts, aes(y=Intensity, x = Sample)) +
    geom_violin( draw_quantiles = 0.5) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
    scale_y_log10() +
    theme_bw() +
    geom_point(data = means, inherit.aes = TRUE)
  
  if(plot){
    if(PDF) pdf(gsub("\\.pdf|$", ".pdf", name), height = 5, width = 7)
    print(p)
    if(PDF) dev.off()
  }else{
    return(p)
  }
}

