# Due to: http://stackoverflow.com/questions/24501245/data-table-throws-object-not-found-error
.datatable.aware=TRUE

#' Filter by sibling peptide correlation
#' @description Filter peptides in traces object based on
#' Sibling Peptide Correlation (SibPepCorr). Since peptides belonging to the same protein
#' should co-elute in a PCP experiment the quality of their signal can be estimated using
#' the average correlation to its sibling peptides. Additionally comparing target and decoy peptides
#' allows for robust FDR estimation.
#' @import data.table
#' @param traces An object of class traces (the trace_type must be "peptide").
#' @param fdr_cutoff Numeric, specifying the FDR cutoff to be applied to select a
#' sibling peptide correlation cutoff
#' @param FFT Numeric, specifying the fraction of false targets (FFT). Default is 1.
#' @param absolute_spcCutoff Numeric, specifying an absolute sibling correlation
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
#' @examples
#' # Load example data
#' tracesRaw <- examplePeptideTraces
#'
#' ## Filter the raw traces object to a protein FDR of 1%
#'   tracesFiltered <- filterBySibPepCorr(traces = tracesRaw,
#'                                        fdr_cutoff = NULL,
#'                                        fdr_type = "protein",
#'                                        FFT = 1,
#'                                        absolute_spcCutoff = 0.2,
#'                                        rm_decoys = FALSE,
#'                                        plot = TRUE,
#'                                        CSV = FALSE)
#'   ## Compare the filtered traces to the raw traces
#' summary(tracesRaw)
#' summary(tracesFiltered)
#' ##--------------------------------------------------------------------------------
#' ## Alternatively the Sibling peptide correlation can also be calculated seperately
#' ##--------------------------------------------------------------------------------
#' # Load example data
#' tracesRaw <- examplePeptideTraces
#'
#' ## Calculate the SibPepCorr of every peptide
#' tracesRawSpc <- calculateSibPepCorr(traces = tracesRaw,
#'                                     plot = TRUE)
#'
#' ## Filter the raw traces object to a protein FDR of 1%
#' tracesFiltered <- filterBySibPepCorr(traces = tracesRawSpc,
#'                                      fdr_cutoff = NULL,
#'                                      fdr_type = "protein",
#'                                      FFT = 1,
#'                                      absolute_spcCutoff = 0.2,
#'                                      rm_decoys = FALSE,
#'                                      plot = TRUE,
#'                                      PDF = FALSE,
#'                                      CSV = FALSE)
#' summary(tracesRaw)
#' summary(tracesFiltered)
#'
filterBySibPepCorr <- function(traces,
                               fdr_cutoff = 0.01,
                               fdr_type = "protein",
                               FFT = 1,
                               absolute_spcCutoff = NULL,
                               rm_decoys = FALSE,
                               plot = TRUE,
                               PDF = FALSE,
                               CSV = FALSE) {

  ## Test traces
  .tracesTest(traces, type = "peptide")
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

  # get spcCutoff from FDR estimation or direct
  decoysContained = length(grep("^DECOY_", traces$trace_annotation$protein_id)) > 0
  if (!is.null(fdr_cutoff)){
    if (decoysContained){
      message("Decoys found...\nEstimating FDR...")
      roctable <- rocSibPepCorr(traces, fdr_type = fdr_type, FFT = FFT, plot = FALSE)
      if(fdr_type == "protein"){
        cutoffForFdr <- roctable[proteinFDR <= fdr_cutoff, min(SibPepCorr_cutoff)]
        fdrReached <- roctable[proteinFDR <= fdr_cutoff, max(proteinFDR)]
        targetProteinsRemaining <- roctable[proteinFDR <= fdr_cutoff, max(n_targetProteins)]
        message("Using SibCorrCutoff of ", round(cutoffForFdr, 4), "\nestimated protein FDR reached: ", round(fdrReached, 4))
      }else if(fdr_type == "peptide"){
        cutoffForFdr <- roctable[peptideFDR <= fdr_cutoff, min(SibPepCorr_cutoff)]
        fdrReached <- roctable[peptideFDR <= fdr_cutoff, max(peptideFDR)]
        targetProteinsRemaining <- roctable[peptideFDR <= fdr_cutoff, max(n_targetProteins)]
        targetPeptidesRemaining <- roctable[peptideFDR <= fdr_cutoff, max(n_targetPeptides)]
        message("Using SibCorrCutoff of ", round(cutoffForFdr, 4), "\nestimated FDR reached: ", round(fdrReached, 4))
      }
      spcCutoff <- cutoffForFdr
    } else if (is.null(absolute_spcCutoff)){
      message("No decoys found...\nPlease supply argument absolute_spcCutoff or label decoy protein_ids with DECOY_ prefix")
    } else{
      message("No Decoys found...Using absulte_spcCutoff: ", absolute_spcCutoff)
      spcCutoff <- absolute_spcCutoff
    }
  } else if (is.null(absolute_spcCutoff)){
    message("No fdr_cutoff or absolute_spcCutoff provided.")
  } else{
    message("No fdr_cutoff provided...Using absulte_spcCutoff: ", absolute_spcCutoff)
    spcCutoff <- absolute_spcCutoff
  }

  ## subset traces object
  tracesFiltered <- traces
  tracesFiltered$trace_annotation <- traces$trace_annotation[SibPepCorr >= spcCutoff]
  tracesFiltered$traces <- traces$traces[traces$trace_annotation$SibPepCorr >= spcCutoff]

  # compare target/decoy numbers before and after
  if (decoysContained){
    ## before
    traces$trace_annotation$DECOY <- 0
    traces$trace_annotation$DECOY[grep('^DECOY_', traces$trace_annotation$protein_id)] <- 1
    #Proteins
    targetProteinsBefore <- traces$trace_annotation[DECOY == 0, length(unique(protein_id))]
    decoyProteinsBefore <- traces$trace_annotation[DECOY == 1, length(unique(protein_id))]
    proteinFdrBefore = decoyProteinsBefore/(targetProteinsBefore) * FFT
    #Peptides
    targetPeptidesBefore <- traces$trace_annotation[DECOY == 0, length(unique(id))]
    decoyPeptidesBefore <- traces$trace_annotation[DECOY == 1, length(unique(id))]
    peptideFdrBefore = decoyPeptidesBefore/(targetPeptidesBefore) * FFT
    ## after
    tracesFiltered$trace_annotation$DECOY <- 0
    tracesFiltered$trace_annotation$DECOY[grep('^DECOY_', tracesFiltered$trace_annotation$protein_id)] <- 1
    #Proteins
    targetProteinsAfter <- tracesFiltered$trace_annotation[DECOY == 0, length(unique(protein_id))]
    decoyProteinsAfter <- tracesFiltered$trace_annotation[DECOY == 1, length(unique(protein_id))]
    protein_fdr_after = decoyProteinsAfter/(targetProteinsAfter) * FFT
    #Peptides
    targetPeptidesAfter <- tracesFiltered$trace_annotation[DECOY == 0, length(unique(id))]
    decoyPeptidesAfter <- tracesFiltered$trace_annotation[DECOY == 1, length(unique(id))]
    peptide_fdr_after = decoyPeptidesAfter/(targetPeptidesAfter) * FFT
  } else{
    # before
    targetProteinsBefore <- traces$trace_annotation[, length(unique(protein_id))]
    decoyProteinsBefore <- 0
    proteinFdrBefore <- NA
    targetPeptidesBefore <- traces$trace_annotation[, length(unique(id))]
    decoyPeptidesBefore <- 0
    peptideFdrBefore <- NA
    # after
    targetProteinsAfter <- tracesFiltered$trace_annotation[, length(unique(protein_id))]
    decoyProteinsAfter <- 0
    protein_fdr_after <- NA
    targetPeptidesAfter <- tracesFiltered$trace_annotation[, length(unique(id))]
    decoyPeptidesAfter <- 0
    peptide_fdr_after <- NA
  }

  ## OUTPUT
  ##-----------------------------------------------------------

  ## Assemble small filter report table
  if(fdr_type == "peptide"){
    SibPepCorrFilterReport <- data.table(cutoff_placeholder = c("pre-filter", "post-filter", "fraction_removed:"),
                                         n_targetProteins = c(targetProteinsBefore, targetProteinsAfter, (1- targetProteinsAfter/targetProteinsBefore)),
                                         n_decoyProteins = c(decoyProteinsBefore, decoyProteinsAfter, (1- decoyProteinsAfter/decoyProteinsBefore)),
                                         protein_fdr = c(proteinFdrBefore, protein_fdr_after, ""),
                                         n_targetPeptides = c(targetPeptidesBefore, targetPeptidesAfter, (1- targetPeptidesAfter/targetPeptidesBefore)),
                                         n_decoy_peptides = c(decoyPeptidesBefore, decoyPeptidesAfter, (1- decoyPeptidesAfter/decoyPeptidesBefore)),
                                         peptide_fdr = c(peptideFdrBefore, peptide_fdr_after, ""))
    setnames(SibPepCorrFilterReport, "cutoff_placeholder", paste0("spcCutoff: ", spcCutoff))
  }else if(fdr_type == "protein"){

    SibPepCorrFilterReport <- data.table(cutoff_placeholder = c("pre-filter", "post-filter", "fraction_removed:"),
                                         n_targetProteins = c(targetProteinsBefore, targetProteinsAfter, (1- targetProteinsAfter/targetProteinsBefore)),
                                         n_decoyProteins = c(decoyProteinsBefore, decoyProteinsAfter, (1- decoyProteinsAfter/decoyProteinsBefore)),
                                         protein_fdr = c(proteinFdrBefore, protein_fdr_after, ""))
    setnames(SibPepCorrFilterReport, "cutoff_placeholder", paste0("spcCutoff: ", spcCutoff))

  }
  ## CSV
  if (CSV){
    if (!is.null(fdr_cutoff)) {
      write.csv(roctable, file = "SibPepCorrFilter_IDFDRtable.csv", row.names = FALSE, quote = FALSE)
    }
    write.csv(SibPepCorrFilterReport, file = "SibPepCorrFilterReport.csv", row.names = FALSE, quote = FALSE)
  }

  # PDF/ROC-like plot
  if (decoysContained & !is.null(fdr_cutoff)){
    if (PDF){
      pdf("SibPepCorrFilter_IDFDRplot.pdf")
    }
    if(plot){
      par(mar=c(5.4,5.4,6.4,5.4))
      if(fdr_type =="protein"){
        xaxis <- roctable$proteinFDR
        targets <- roctable$n_targetProteins
        true_targets <- roctable$n_true_targetProteins
        targets_remaining <- targetProteinsRemaining
      }else if(fdr_type =="peptide"){
        xaxis <- roctable$peptideFDR
        targets <- roctable$n_targetPeptides
        true_targets <- roctable$n_true_targetPeptides
        targets_remaining <- targetPeptidesRemaining
      }

      plot(xaxis, targets, type = "l", lty = 2,
           lwd = 2,
           main = paste0("ID-FDR Curve\n spcCutoff: ", round(spcCutoff, 4), "\ntarget proteins remaining: ",
                         targetProteinsRemaining, "\nestimated FDR: ", round(fdrReached, 4)),
           xlim = c(0,0.1), ylim = c(0, 1.02* max(targets)), xlab='FDR',ylab='n')
      lines(xaxis, true_targets, lty = 1, lwd = 2)
      abline(h = targets_remaining, lty = 2, col = "red")
      abline(v = fdrReached, lty = 2, col = "red")
      points(fdrReached, targets_remaining, pty = 2, lwd = 3, col = "red")
      legend("bottomright", lty = c(2,1), legend = c("all target proteins", "true target proteins"))
    }
    if (PDF){
      dev.off()
    }
  }

  # remove decoys
  if (rm_decoys) {
    n_decoys <- length(grep("^DECOY", tracesFiltered$trace_annotation$protein_id))
    if ( n_decoys > 0){
      idx_decoys <- grep("^DECOY_",tracesFiltered$trace_annotation$protein_id)
      tracesFiltered$traces <- tracesFiltered$traces[-idx_decoys]
      tracesFiltered$trace_annotation<- tracesFiltered$trace_annotation[-idx_decoys]
      message(n_decoys, " decoys removed")
    } else {
      message("no decoys contained/removed")
    }
  }

  #return the filtered traces data
  message("Proteins remaining in dataset: ", targetProteinsAfter)
  return(tracesFiltered)

}
