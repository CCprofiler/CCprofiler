#' Calculates the correlation of every petide between replicates
#' @description Average pairwise pearson correlation of the intensity vectors between
#' replicates.
#' @param traces Object of class \code{tracesList}.
#' @param design_matrix data.table, describing the \code{tracesList} object.
#' @param compare_within Character string, name of the column defining the conditions.
#' @param add Logical,add repPepCorr to \code{tracesList} object or return a table.
#' @return Object of class traces with additional repPepCorr column or data.table
#' contaioning repPepCors.
#' @export

calculateRepPepCorr <- function(traces, design_matrix, compare_within = NULL,
                                add = T){
  
  # Check input type
  
  if (any(sapply(traces, "[[", "trace_type") != "peptide")){
    stop("Replicate peptide correlation can only be calculated on traces of type peptide")
  }
  .tracesListTest(traces, type = "peptide")
  for(i in 1:length(traces)){
    if(any(names(traces[[i]]$trace_annotation)== "RepPepCorr" & add == T)){
      traces[[i]]$trace_annotation[,RepPepCorr := NULL]
      message(paste0("Found RepPepCorr in ", names(traces)[i], ": Overwriting..."))
    }
    setkey(traces[[i]]$traces, "id")
  }
  
  if(!is.null(compare_within)){
    compare_col <- which(names(design_matrix) == compare_within)
    stopifnot(length(compare_col) == 1)
    conditions <- unique(design_matrix[[compare_within]])
  }else{
    design_matrix$Condition <- "all Conditions"
    conditions <- "all Conditions"
    compare_within <- "Condition"
  }
  
  pwcorr <- calculatePairwiseRepPepCorr(traces)
  
  if(add){
    for(condition in conditions){
      samplesInCond <- design_matrix[get(compare_within) == condition, Sample_name]
      message(paste0("Calculating RepPepCorrs for ", condition, "..."))
      pc <- pwcorr[Sample1 %in% samplesInCond & Sample2 %in% samplesInCond & !is.na(RepPepCorr)]
      # pc[, detected_in := .N, by = id]
      for(sample in samplesInCond){
        pcs <- pc[(Sample1 == sample | Sample2 == sample)]
        meanpc <- pcs[, .(RepPepCorr = mean(RepPepCorr)), by = id]
        traces[[sample]]$trace_annotation <- merge(traces[[sample]]$trace_annotation,
                                                   meanpc,
                                                   by = "id", all.x = T)
      }
    }
    .tracesListTest(traces)
    return(traces)
  }else{
    return(pwcorr)      
  }
}

# calculateRepPepCorr <- function(traces, design_matrix, compare_within = NULL,
#                                 add = T){
#   
#   # Check input type
#   
#   if (any(sapply(traces, "[[", "trace_type") != "peptide")){
#     stop("Replicate peptide correlation can only be calculated on traces of type peptide")
#   }
#   .tracelListTest(traces, type = "peptide")
#   for(i in 1:length(traces)){
#     if(any(names(traces[[i]]$trace_annotation)== "RepPepCorr" & add == T)){
#       traces[[i]]$trace_annotation[,RepPepCorr := NULL]
#       message(paste0("Found RepPepCorr in ", names(traces)[i], ": Overwriting..."))
#     }
#     setkey(traces[[i]]$traces, "id")
#   }
#   
#   if(!is.null(compare_within)){
#     compare_col <- which(names(design_matrix) == compare_within)
#     stopifnot(length(compare_col) == 1)
#     conditions <- unique(design_matrix[[compare_within]])
#   }else{
#     design_matrix$Condition <- "all Conditions"
#     conditions <- "all Conditions"
#     compare_within <- "Condition"
#   }
#   
#   repPepCorr <- lapply(conditions, function(condition){
#     sample_names <- design_matrix[get(compare_within) == condition, Sample_name]
#     
#     message(paste0("Calculating RepPepCorr for Replicates within ", condition, "..."))
#     df <- do.call("rbind", lapply(traces[sample_names], function(x) x$traces))
#     setkey(df, "id")
#     # Subset to intersection peptides
#     # if(is.null(intersect)){
#     #   # intersect_peps <- Reduce(intersect, lapply(traces[sample_names], function(x) x$trace_annotation$id))
#     #   # df <- df[intersect]
#     #   intersect <- length(sample_names)
#     # }
#     # df <- df[,if(.N >=intersect) .SD,by=id]
#     # Calculate the mean correllation of the intersection peptides
#     df_cors <- df[,
#                   .({c <- cor(t(.SD))
#                   mean(c[lower.tri(c)])},
#                   nrow(c)),
#                   by = id]
#     names(df_cors) <- c("id", "RepPepCorr", "detected_in")
#     # sapply(intersect, function(pep){
#     #   df_p <- df
#     #   df[, id:= NULL]
#     #   df_cor <- cor(t(df))
#     #   mean(df_cor[lower.tri(df_cor)])
#     # }) 
#     # 
#     df_cors
#   })
#   names(repPepCorr) <- conditions
#   
#   if(add){
#     for(condition in conditions){
#       sample_names <- design_matrix[design_matrix[[2]] == condition, Sample_name]
#       for(sample_name in sample_names){
#         traces[[sample_name]]$trace_annotation <- merge(traces[[sample_name]]$trace_annotation,
#                                                         repPepCorr[[condition]],
#                                                         by = "id", all.x = T)
#       }
#     }
#     return(traces)
#   }else{
#     return(repPepCorr)      
#   }
# }

calculatePairwiseRepPepCorr <- function(traces, comparisons = NULL){
  
  # Check input type
  if (any(sapply(traces, "[[", "trace_type") != "peptide")){
    stop("Replicate peptide correlation can only be calculated on traces of type peptide")
  }
  # for(i in 1:length(traces)){
  #   setkey(traces[[i]]$traces, "id")
  # }
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

#' plot RepPepCorrDensities
#' @description Plot replicate peptide correlation in traces object of type peptide.
#' @import data.table
#' @param traces An object of type traces.
#' @param PDF logical, wether to print RepPepCorr density plot to a PDF file. Deafaults to \code{FALSE}.
#' @return Plot.
#' @export

plotRepPepCorrDensities <- function(traces, PDF = FALSE,
                                    name = "ReplicatePeptideCorrelationDensity"){
  UseMethod("plotRepPepCorrDensities", traces)
}

#' @describeIn plotRepPepCorrDensities
plotRepPepCorrDensities.traces <- function(traces, PDF = FALSE,
                                           name = "ReplicatePeptideCorrelationDensity"){
  
  ## Test traces
  .tracesTest(traces, type = "peptide")
  
  # check whether decoys are present in the input peptide traces object
  trace_annotation <- traces$trace_annotation
  decoys_present = TRUE
  if (sum(grep("^DECOY_", trace_annotation$protein_id)) == 0){
    message("No decoy entries in trace_annotation$protein_id column \n
            No decoy density will be plotted")
    decoys_present = FALSE
  }
  
  #check whether RepPepCorr has been calculated/is contained in trace_annotation
  if (!("RepPepCorr" %in% names(traces$trace_annotation))){
    stop("No RepPepCorr has been calculated on this dataset. Use calculateRepPepCorr function.")
  }
  
  dens_all <- density(na.omit(trace_annotation$RepPepCorr))
  dens_targets <- density(na.omit(trace_annotation[grep("^DECOY", trace_annotation$protein_id,invert = TRUE),]$RepPepCorr))
  
  if(PDF == TRUE){
    pdf(gsub("\\.pdf$|$", ".pdf", name))
  }
  
  plot(dens_targets$x, dens_targets$y*dens_targets$n, lty = 1, lwd = 3,
       type = "l", ylab = "scaled frequency", xlab = "RepPepCorr",
       main = name)
  
  if (decoys_present){
    dens_decoys <- density(na.omit(trace_annotation[grep("^DECOY_", trace_annotation$protein_id),]$RepPepCorr))
    lines(dens_decoys$x, dens_decoys$y*dens_decoys$n, lty = 2, lwd = 3)
    legend("topleft", lty = c(1,3), lwd = c(3,3), legend = c("target peptides", "decoy peptides"))
  } else{
    legend("topleft", lty = c(1), lwd = c(3), legend = c("target peptides"))
  }
  
  if(PDF == TRUE){
    dev.off()
  }
}

#' @describeIn plotRepPepCorrDensities
plotRepPepCorrDensities.tracesList <- function(traces, PDF = FALSE,
                                               name = "ReplicatePeptideCorrelationDensity"){
  .tracesListTest(traces)
  
  if(PDF) pdf(gsub("$|\\.pdf$", ".pdf",name))
  res <- lapply(names(traces), function(tr){
    plotRepPepCorrDensities.traces(traces = traces[[tr]],
                                   PDF = F,
                                   name = paste(name, tr, sep = " - "))
  })
  if(PDF) dev.off()
}