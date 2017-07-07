# Due to: http://stackoverflow.com/questions/24501245/data-table-throws-object-not-found-error
.datatable.aware=TRUE

#' plot SibPepCorrDensities
#' @description Plot sibling peptide correlation in traces object of type peptide.
#' @import data.table
#' @param traces An object of type traces.
#' @param PDF logical, wether to print SibPepCorr density plot to a PDF file. Deafaults to \code{FALSE}.
#' @return Plot.
#' @export
#' @example 
#' ## Load example data
#' peptideTraces <- examplePeptideTracesFiltered
#' 
#' ## Plot the SibPepCorr density plot
#' plotSibPepCorrDensities(peptideTraces)

plotSibPepCorrDensities <- function(traces, PDF = FALSE){
  
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

  #check whether SibPepCorr has been calculated/is contained in trace_annotation
  if (!("SibPepCorr" %in% names(traces$trace_annotation))){
    stop("No SibPepCorr has been calculated on this dataset. Use calculateSibPepCorr function.")
  }

  dens_all <- density(na.omit(trace_annotation$SibPepCorr))
  dens_targets <- density(na.omit(trace_annotation[grep("^DECOY", trace_annotation$protein_id,invert = TRUE),]$SibPepCorr))

  if(PDF == TRUE){
    pdf("CorrFilter_SibPepCorr_distributions_target_decoy.pdf")
  }

  plot(dens_targets$x, dens_targets$y*dens_targets$n, lty = 1, lwd = 3,
       type = "l", ylab = "scaled frequency", xlab = "SibPepCorr",
       main = "Sibling Peptide Correlation Density")

  if (decoys_present){
    dens_decoys <- density(na.omit(trace_annotation[grep("^DECOY_", trace_annotation$protein_id),]$SibPepCorr))
    lines(dens_decoys$x, dens_decoys$y*dens_decoys$n, lty = 2, lwd = 3)
    legend("topleft", lty = c(1,3), lwd = c(3,3), legend = c("target peptides", "decoy peptides"))
  } else{
    legend("topleft", lty = c(1), lwd = c(3), legend = c("target peptides"))
  }

  if(PDF == TRUE){
    dev.off()
  }
}
