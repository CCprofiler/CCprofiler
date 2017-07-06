# Due to: http://stackoverflow.com/questions/24501245/data-table-throws-object-not-found-error
.datatable.aware=TRUE

#' Filter by sibling peptide correlation
#' @description Filter peptides in traces object based on
#' Sibling Peptide Correlation (SibPepCorr). Since peptides belonging to the same protein
#' should co-elute in a PCP experiment the quality of their signal can be estimated using
#' the average correlation to its sibling peptides. Additionally comparing target and decoy peptides
#' allows for robust FDR estimation.
#' @details 
#' @import data.table
#' @param traces An object of class traces (the trace_type must be "peptide").
#' @param fdr_cutoff Numeric, specifying the FDR cutoff to be applied to select a
#' sibling peptide correlation cutoff
#' @param FFT Numeric, specifying the fraction of false targets (FFT). Default is 1.
#' @param absolute_spc_cutoff Numeric, specifying an absolute sibling correlation
#' cutoff to be applied. Only used if fdr_cutoff is not provided. Default is NULL.
#' @param rm_decoys Logical, specifying if decoys should be removed.
#' It is strongly recommended to keep the decoys, for visual assessment. Default is FALSE.
#' @param plot Logical, wether to print Score distribution and precision-recall plots
#' to R console.
#' @param PDF Logical, wether to print Score distribution and precision-recall plots
#' to "SibPepCorrFilter_IDFDRplot.pdf" in working directory
#' @param CSV Logical, write a table output of the FDR estimation results to the working directory
#' @return a traces object, only containing peptides passing the specified spc/fdr cutoff. If the 
#' input traces object does not contain calculated spc values for every peptide a SibPepCorr column 
#' is added.
#' @export
#' 
## Load example data
# tracesRaw <- examplePeptideTraces
# 
# ## Filter the raw traces object to a protein FDR of 1%
# tracesFiltered <- filterBySibPepCorr(traces = tracesRaw,
#                                      fdr_cutoff = 0.01,
#                                      fdr_type = "peptide",
#                                      FFT = 0.4,
#                                      absolute_spc_cutoff = NULL,
#                                      rm_decoys = FALSE,
#                                      plot = TRUE)
# ## Compare the filtered traces to the raw traces
# summary(tracesRaw)
# summary(tracesFiltered)
# ##--------------------------------------------------------------------------------
# ## Alternatively the Sibling peptide correlation can also be calculated seperately
# ##--------------------------------------------------------------------------------
# tracesRaw <- examplePeptideTraces
# 
# tracesRawSpc <- calculateSibPepCorr(traces = tracesRaw,
#                                     plot = TRUE)
# 
# tracesFiltered <- filterBySibPepCorr(traces = tracesRawSpc,
#                                      fdr_cutoff = 0.01,
#                                      fdr_type = "protein",
#                                      FFT = 0.4,
#                                      absolute_spc_cutoff = NULL,
#                                      rm_decoys = FALSE,
#                                      plot = TRUE,
#                                      PDF = FALSE,
#                                      CSV = FALSE)
# summary(tracesRaw)
# summary(tracesFiltered)
# 
filterBySibPepCorr <- function(traces,
                               fdr_cutoff = 0.01,
                               fdr_type = "protein",
                               FFT = 1,
                               absolute_spc_cutoff = NULL,
                               rm_decoys = FALSE,
                               plot = TRUE,
                               PDF = FALSE,
                               CSV = FALSE) {
  
  ## If sibling peptide correlation hass not been calculated yet, do it
  if (!("SibPepCorr" %in% names(traces$trace_annotation))){
    message("Sibling peptide correlation not yet calculated for this dataset\nCalculating SibPepCorr(spc)...")
    traces <- calculateSibPepCorr(traces, plot = plot, PDF = PDF)
  } else{
    message("Sibling peptide correlation values found...")
  }
  if(fdr_type == "peptide"){
    message("Warning: It is strongly recommended to filter with protein FDR when using the complex
            centric workflow. Please proceed with caution.")
  }else if(fdr_type != "protein") stop("Parameter fdr_type must be 'protein' or 'peptide'")
  
  # get spc_cutoff from FDR estimation or direct
  decoys_contained = length(grep("^DECOY_", traces$trace_annotation$protein_id)) > 0
  if (decoys_contained){
    message("Decoys found...\nEstimating FDR...")
    roctable <- rocSibPepCorr(traces, fdr_type = fdr_type, FFT = FFT, plot = FALSE)
    if(fdr_type == "protein"){
      cutoff_for_fdr <- roctable[proteinFDR <= fdr_cutoff, min(SibPepCorr_cutoff)]
      fdr_reached <- roctable[proteinFDR <= fdr_cutoff, max(proteinFDR)]
      target_proteins_remaining <- roctable[proteinFDR <= fdr_cutoff, max(n_target_proteins)]
      message("Using SibCorrCutoff of ", round(cutoff_for_fdr, 4), "\nestimated protein FDR reached: ", round(fdr_reached, 4))
    }else if(fdr_type == "peptide"){
      cutoff_for_fdr <- roctable[peptideFDR <= fdr_cutoff, min(SibPepCorr_cutoff)]
      fdr_reached <- roctable[peptideFDR <= fdr_cutoff, max(peptideFDR)]
      target_proteins_remaining <- roctable[peptideFDR <= fdr_cutoff, max(n_target_proteins)]
      target_peptides_remaining <- roctable[peptideFDR <= fdr_cutoff, max(n_target_peptides)]
      message("Using SibCorrCutoff of ", round(cutoff_for_fdr, 4), "\nestimated FDR reached: ", round(fdr_reached, 4))
    }
    spc_cutoff <- cutoff_for_fdr
  } else if (is.null(absolute_spc_cutoff)){
    message("No decoys found...\nPlease supply argument absolute_spc_cutoff or label decoy protein_ids with DECOY_ prefix")
  } else{
    message("No Decoys found...Using absulte_spc_cutoff: ", absolute_spc_cutoff)
    spc_cutoff <- absolute_spc_cutoff
  }
  
  # subset traces object
  traces.filtered <- traces
  traces.filtered$trace_annotation <- traces$trace_annotation[SibPepCorr >= spc_cutoff]
  traces.filtered$traces <- traces$traces[traces$trace_annotation$SibPepCorr >= spc_cutoff]
  
  # compare target/decoy numbers before and after
  if (decoys_contained){
    # before
    traces$trace_annotation$DECOY <- 0
    traces$trace_annotation$DECOY[grep('^DECOY_', traces$trace_annotation$protein_id)] <- 1
    target_proteins_before <- traces$trace_annotation[DECOY == 0, length(unique(protein_id))]
    decoy_proteins_before <- traces$trace_annotation[DECOY == 1, length(unique(protein_id))]
    protein_fdr_before = decoy_proteins_before/(decoy_proteins_before+target_proteins_before) * FFT
    # after
    traces.filtered$trace_annotation$DECOY <- 0
    traces.filtered$trace_annotation$DECOY[grep('^DECOY_', traces.filtered$trace_annotation$protein_id)] <- 1
    target_proteins_after <- traces.filtered$trace_annotation[DECOY == 0, length(unique(protein_id))]
    decoy_proteins_after <- traces.filtered$trace_annotation[DECOY == 1, length(unique(protein_id))]
    protein_fdr_after = decoy_proteins_after/(decoy_proteins_after+target_proteins_after) * FFT
  } else{
    # before
    target_proteins_before <- traces$trace_annotation[, length(unique(protein_id))]
    decoy_proteins_before <- 0
    protein_fdr_before <- NA
    # after
    target_proteins_after <- traces.filtered$trace_annotation[, length(unique(protein_id))]
    decoy_proteins_after <- 0
    protein_fdr_after <- NA
  }
  
  # OUTPUT
  # Assemble small filter report table
  SibPepCorrFilter_report <- data.table(cutoff_placeholder = c("pre-filter", "post-filter", "fraction_removed:"),
                                        n_target_proteins = c(target_proteins_before, target_proteins_after, (1- target_proteins_after/target_proteins_before)),
                                        n_decoy_proteins = c(decoy_proteins_before, decoy_proteins_after, (1- decoy_proteins_after/decoy_proteins_before)),
                                        protein_fdr = c(protein_fdr_before, protein_fdr_after, ""))
  setnames(SibPepCorrFilter_report, "cutoff_placeholder", paste0("spc_cutoff: ", spc_cutoff))
  
  # CSV
  if (CSV){
    write.csv(roctable, file = "SibPepCorrFilter_IDFDRtable.csv", row.names = FALSE, quote = FALSE)
    write.csv(SibPepCorrFilter_report, file = "SibPepCorrFilter_report.csv", row.names = FALSE, quote = FALSE)
  }
  
  # PDF/ROC-like plot
  if (decoys_contained){
    if (PDF){
      pdf("SibPepCorrFilter_IDFDRplot.pdf")
    }
    if(plot){
      par(mar=c(5.4,5.4,6.4,5.4))
      if(fdr_type =="protein"){
        xaxis <- roctable$proteinFDR
        targets <- roctable$n_target_proteins
        true_targets <- roctable$n_true_target_proteins
        targets_remaining <- target_proteins_remaining
      }else if(fdr_type =="peptide"){
        xaxis <- roctable$peptideFDR
        targets <- roctable$n_target_peptides
        true_targets <- roctable$n_true_target_peptides
        targets_remaining <- target_peptides_remaining
      }
      
      plot(xaxis, targets, type = "l", lty = 2,
           lwd = 2,
           main = paste0("ID-FDR Curve\n spc_cutoff: ", round(spc_cutoff, 4), "\ntarget proteins remaining: ",
                         target_proteins_remaining, "\nestimated FDR: ", round(fdr_reached, 4)),
           xlim = c(0,0.1), ylim = c(0, 1.02* max(targets)), xlab='FDR',ylab='n')
      lines(xaxis, true_targets, lty = 1, lwd = 2)
      abline(h = targets_remaining, lty = 2, col = "red")
      abline(v = fdr_reached, lty = 2, col = "red")
      points(fdr_reached, targets_remaining, pty = 2, lwd = 3, col = "red")
      legend("bottomright", lty = c(2,1), legend = c("all target proteins", "true target proteins"))
    }
    if (PDF){
      dev.off()
    }
  }
  
  # remove decoys
  if (rm_decoys) {
    n_decoys <- length(grep("^DECOY", traces.filtered$trace_annotation$protein_id))
    if ( n_decoys > 0){
      idx_decoys <- grep("^DECOY_",traces.filtered$trace_annotation$protein_id)
      traces.filtered$traces <- traces.filtered$traces[-idx_decoys]
      traces.filtered$trace_annotation<- traces.filtered$trace_annotation[-idx_decoys]
      message(n_decoys, " decoys removed")
    } else {
      message("no decoys contained/removed")
    }
  }
  
  #return the filtered traces data
  message("Proteins remaining in dataset: ", target_proteins_after)
  return(traces.filtered)
  
}
