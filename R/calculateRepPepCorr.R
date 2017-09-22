

calculateRepPepCorr <- function(traces, design_matrix, compare_within = NULL,
                                add = T, sa = NULL){
  
  # Check input type
  if (any(sapply(traces, "[[", "trace_type") != "peptide")){
    stop("Replicate peptide correlation can only be calculated on traces of type peptide")
  }
  for(i in 1:length(traces)){
    if(any(names(traces[[i]]$trace_annotation)== "RepPepCorr" & add == T)){
      traces[[i]]$trace_annotation[,RepPepCorr := NULL]
      message(paste0("Found RepPepCorr in ", names(traces)[i], ": Overwriting..."))
    }
    setkey(traces[[i]]$traces, "id")
  }
  
  
  # Find peptides present in all replicates of the conditions
  if(!is.null(compare_within)){
    compare_col <- which(names(design_matrix) == compare_within)
    stopifnot(length(compare_col) == 1)
    conditions <- unique(design_matrix[[compare_within]])
  }else{
    design_matrix$Condition <- "all Conditions"
    conditions <- "all Conditions"
    compare_within <- "Condition"
  }
  
  repPepCorr <- lapply(conditions, function(condition){
    sample_names <- design_matrix[get(compare_within) == condition, Sample_name]
    
    message(paste0("Calculating RepPepCorr for Replicates within ", condition, "..."))
    df <- do.call("rbind", lapply(traces[sample_names], function(x) x$traces))
    setkey(df, "id")
    # Subset to intersection peptides
    # if(is.null(intersect)){
    #   # intersect_peps <- Reduce(intersect, lapply(traces[sample_names], function(x) x$trace_annotation$id))
    #   # df <- df[intersect]
    #   intersect <- length(sample_names)
    # }
    # df <- df[,if(.N >=intersect) .SD,by=id]
    # Calculate the mean correllation of the intersection peptides
    df_cors <- df[,
                  .({c <- cor(t(.SD))
                  mean(c[lower.tri(c)])},
                  nrow(c)),
                  by = id]
    names(df_cors) <- c("id", "RepPepCorr", "detected_in")
    # sapply(intersect, function(pep){
    #   df_p <- df
    #   df[, id:= NULL]
    #   df_cor <- cor(t(df))
    #   mean(df_cor[lower.tri(df_cor)])
    # }) 
    # 
    df_cors
  })
  names(repPepCorr) <- conditions
  
  if(add){
    for(condition in conditions){
      sample_names <- design_matrix[design_matrix[[2]] == condition, Sample_name]
      for(sample_name in sample_names){
        traces[[sample_name]]$trace_annotation <- merge(traces[[sample_name]]$trace_annotation,
                                                        repPepCorr[[condition]],
                                                        by = "id", all.x = T)
      }
    }
    return(traces)
  }else{
    return(repPepCorr)      
  }
}

