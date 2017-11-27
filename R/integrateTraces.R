#' Add or average traces objects to generate consensus traces.
#' @param traces object of class tracesList.
#' (must be a valid column of the design matrix). If \code{NULL}, everything is collapsed.
#' @param min_present Numeric, the minimum number of samples in which a trace must be detected.
#' @param aggr_fun Character string, the method of aggregating the traces.
#'  Must be one of c("sum", "mean").
#' @return tracesList object containing the integrated traces
#' @export

integrateTraceIntensities <- function(traces,
                                      design_matrix = NULL,
                                      integrate_within = NULL,
                                      min_present = 1, 
                                      aggr_fun = c("sum", "mean")){
  
  for(i in 1:length(traces)){
    .tracesTest(traces[[i]], type = traces[[1]]$trace_type)
  }
  aggr_fun <- match.arg(aggr_fun)
  
  if(!is.null(integrate_within)){
    if(is.null(design_matrix)){
      stop("If integrate_within is specified, you have to specify a design matrix as well")
    }
    compare_col <- which(names(design_matrix) == integrate_within)
    stopifnot(length(compare_col) == 1)
    conditions <- unique(design_matrix[[integrate_within]])
  }else{
    design_matrix <- data.table(Sample_name = names(traces),
                                Condition = "all Conditions")
    conditions <- "all Conditions"
    integrate_within <- "Condition"
  }
  oldTraces <- traces
  tracesRes <- list()
  for(condition in conditions){
    samples <- design_matrix[get(integrate_within) == condition, Sample_name]
    traces <- oldTraces[samples]
    
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
    tracesRes[[condition]] <- resTraces
  }
  
  if(length(tracesRes) >1){
    class(tracesRes) <- "tracesList"
    .tracesListTest(tracesRes)
  }else{
    tracesRes <- tracesRes[[1]]
    .tracesTest(tracesRes)
  }
  return(tracesRes)
}
