#' annotateMassDistribution
#' @description Estimates fraction of assembled vs. monomeric mass per id.
#' @param traces An object of type traces. These should be the raw, unprocessed traces.
#' @export
annotateMassDistribution <- function(traces){
  UseMethod("annotateMassDistribution", traces)
}

#' @describeIn annotateMassDistribution Estimates fraction of assembled vs. monomeric mass per id.
#' @export
annotateMassDistribution.traces <- function(traces){
  .tracesTest(traces)

  if(! ("protein_mw" %in% names(traces$trace_annotation))) {
    stop("annotateMassDistribution requires information on the monomer molecular weight. Please see annotateTraces.")
  }

  if(! ("molecular_weight" %in% names(traces$fraction_annotation))) {
    stop("annotateMassDistribution requires information on the molecular weight assigned to each fraction. Please see annotateMolecularWeight.")
  }

  intMat <- getIntensityMatrix(traces)
  mw_info <- traces$trace_annotation
  max_trace <- max(traces$fraction_annotation$id)
  mw_info[,assembly_boundary := 2*protein_mw]
  mw_info[is.na(assembly_boundary), assembly_boundary := 0]
  mw_info[,assembly_boundary_fraction := unlist(lapply(assembly_boundary,function(x){max(traces$fraction_annotation[molecular_weight>=x]$id)}))]

  # Get intensity sum
  mw_info[,sum_assembled := sum(intMat[id,1:assembly_boundary_fraction-1]), id]
  mw_info[,sum_monomeric := sum(intMat[id,assembly_boundary_fraction:max_trace]), id]
  mw_info[,sum := sum(intMat[id,]), id]
  mw_info[is.na(assembly_boundary), sum_assembled := 0]
  mw_info[is.na(assembly_boundary), sum_monomeric := 0]
  mw_info[, sum_assembled_norm := 1/sum*sum_assembled]
  mw_info[, sum_monomeric_norm := 1/sum*sum_monomeric]
  #mw_info[, test := sum_monomeric_norm+sum_assembled_norm]

  traces$trace_annotation <- mw_info
  .tracesTest(traces)
  return(traces)
}

#' @describeIn annotateMassDistribution Estimates fraction of assembled vs. monomeric mass per id.
#' @export
annotateMassDistribution.tracesList <- function(tracesList){
  .tracesListTest(tracesList)
  tracesListRes <- lapply(tracesList, annotateMassDistribution.traces)
  class(tracesListRes) <- "tracesList"
  .tracesListTest(tracesListRes)
  return(tracesListRes)
}

#' getMassAssemblyChange
#' @description Estimates change of assembled vs. monomeric mass per id.
#' @param tracesList An object of type traces.list including assembled mass 
#' annotation as produced by \code{annotateMassDistribution}.
#' @param design_matrix data.table, design matrix describing the architecture of the tracesList object.
#' @export
getMassAssemblyChange <- function(tracesList, design_matrix){
  .tracesListTest(tracesList)
  samples <- unique(design_matrix$Sample)
  if(length(samples) != 2) {
    stop("This function is only available for comparing 2 conditions.")
  }
  if(! all(samples %in% names(tracesList))) {
    stop("tracesList and design_matrix do not match. Pleas check sample names.")
  }
  if (! "sum_assembled_norm" %in% names(tracesList[[tr]]$trace_annotation)) {
    stop("No assembled mass annotation available, please run annotateMassDistribution first.")
  }
  res <- lapply(names(tracesList), function(tr){
    vals <- subset(tracesList[[tr]]$trace_annotation,select=c("protein_id","sum_assembled_norm"))
    vals[,Sample := tr]
    return(vals)
  })
  res <- do.call(rbind, res)
  res_cast <- dcast(res, formula = protein_id ~ Sample, value.var=c("sum_assembled_norm"))
  res_cast[, change := log2(get(samples[1])/(get(samples[2])))]
  res_cast[is.nan(change), change := 0]
  res_cast[, testOrder := paste0(samples[1],".vs.",samples[2])]
  return(res_cast[])
}

#' plotMassAssemblyChange between 2 conditions
#' @param assamblyTest data.table, a data.table with test statistics.
#' An assamblyTest can be produced with \code{getMassAssemblyChange}.
#' @param FC_cutoff Numeric fold change cutoff, default is 2.
#' @param name character string specifying the name of output if PDF=TRUE, default is "massAssemblyChange".
#' @param PDF logical if PDF should be created, default is FALSE.
#' @return plot
#' @export
plotMassAssemblyChange <- function(assamblyTest, FC_cutoff=2, name="massAssemblyChange", PDF=FALSE){
  if (PDF) {
    pdf(paste0(name,".pdf"))
  }
  no_change <- nrow(subset(assamblyTest, abs(change) <= FC_cutoff))
  more_assembled <- nrow(subset(assamblyTest, change < -FC_cutoff))
  less_assembled <- nrow(subset(assamblyTest, change > FC_cutoff))
  pie(c(no_change,more_assembled,less_assembled),
      labels=c(paste("abs(logFC) < ",FC_cutoff,"\n",no_change),
               paste0("logFC < -",FC_cutoff,"\n",more_assembled),
               paste0("logFC > ",FC_cutoff,"\n",less_assembled)),
      main = paste0(unique(assamblyTest$testOrder)))
  
  h <- ggplot(assamblyTest, aes(x=change)) + geom_histogram(binwidth = 1) + theme_classic()
  print(h)
  if (PDF) {
    dev.off()
  }
}

