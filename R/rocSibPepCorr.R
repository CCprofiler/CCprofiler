# Due to: http://stackoverflow.com/questions/24501245/data-table-throws-object-not-found-error
.datatable.aware=TRUE

#' ROC.SibPepCorr
#' @description Calculate ROC characteristics based on sibling peptide correlation.
#' @import data.table
#' @param traces An object of type \code{traces.obj}.
#' @param FFT Numeric specifying the fraction of false targets (FFT). Default is 1.
#' @param stepsize Numeric specifying the number of steps for calculating ROC characteristics. Default is 0.001.
#' @param plot Logical, whether to generate a plot. Default is \code{TRUE}.
#' @param PDF Logical, whether to save plot as PDF. Default is \code{FALSE}.
#' @param CSV Logical, whether to save ROC characteristics in a CSV file. Default is \code{FALSE}.
#' @param name Character string with the name of the saved PDF or CSV file.
#' Only used if PDF or CSV is set to \code{TRUE}. Default is "SibPepCorr_ROC".
#' @return A data.table containing ROC characteristics.
#' @export
#' @examples
#' # Load example data
#' tracesRaw <- examplePeptideTraces
#' ## Calculate the sibPepCorr
#' tracesRawSpc <- calculateSibPepCorr(traces = tracesRaw,
#'                                       plot = TRUE)
#' ## Calculate protein FDR at different correllation cutoffs
#' fdrTable <- rocSibPepCorr(traces = tracesRawSpc,
#'                           plot = TRUE,
#'                           fdr_type = "protein")
#'
#' head(fdrTable)

rocSibPepCorr <- function(traces,
                          fdr_type = "protein",
                          FFT = 1,
                          stepsize = 0.001,
                          plot = TRUE,
                          PDF = FALSE,
                          CSV = FALSE,
                          name = "SibPepCorr_ROC"){
  ## Test traces
  .tracesTest(traces, type = "peptide")

  ## check whether SibPepCorr has been calculated/is contained in trace_annotation
  if (!("SibPepCorr" %in% names(traces$trace_annotation))){
    stop("No SibPepCorr has been calculated on this dataset. Use calculateSibPepCorr function.")
  }
  ##  Make sure that teh right FDR type is specified
  if(!(fdr_type %in% c("protein", "peptide"))) stop("Parameter fdr_type must be 'protein' or 'peptide'")

  trace_annotation <- traces$trace_annotation
  trace_annotation$DECOY <- 0
  trace_annotation$DECOY[grep('^DECOY_', trace_annotation$protein_id)] <- 1
  corr_range <- sort(unique(trace_annotation$SibPepCorr))
  corrTest <- seq(min(corr_range), max(corr_range), stepsize)
  ncorrTest <- length(corrTest)
  targetProteins <- numeric(length = ncorrTest)
  decoyProteins <- numeric(length = ncorrTest)
  targetPeptides <- numeric(length = ncorrTest)
  decoyPeptides <- numeric(length = ncorrTest)

  # Step through the spc cutoffs and calculate the FDR

  if(fdr_type == "peptide"){
    for (i in 1:ncorrTest) {
      remainingTraces <- trace_annotation[SibPepCorr >= corrTest[i]]
      ## Count proteins
      targetProts <- unique(remainingTraces[DECOY == 0]$protein_id)
      targetProteins[i] <- length(targetProts)
      decoyProts <- unique(remainingTraces[DECOY == 1]$protein_id)
      decoyProteins[i] <- length(decoyProts)

      ## Count peptides
      targetPeptides[i] <- nrow(remainingTraces[DECOY == 0])
      decoyPeptides[i] <- nrow(remainingTraces[DECOY == 1])
    }
  }else if(fdr_type == "protein"){
    for (i in 1:ncorrTest) {
      remainingTraces <- trace_annotation[SibPepCorr >= corrTest[i]]
      ## Count proteins
      targetProts <- unique(remainingTraces[DECOY == 0]$protein_id)
      targetProteins[i] <- length(targetProts)
      decoyProts <- unique(remainingTraces[DECOY == 1]$protein_id)
      decoyProteins[i] <- length(decoyProts)
    }
  }
## Estimate the FDR
  if(fdr_type == "peptide"){

  fdr_protein <- pmin(1, FFT*decoyProteins/(targetProteins))
  fdr_peptides <- pmin(1, FFT*decoyPeptides/(targetPeptides))
  true_targetProteins <- as.integer(targetProteins * (1-fdr_protein))
  true_targetPeptides <- as.integer(targetPeptides * (1-fdr_peptides))
  resulttable <- as.data.table(cbind(corrTest,
                                     targetProteins, true_targetProteins, decoyProteins,fdr_protein,
                                     targetPeptides, true_targetPeptides, decoyPeptides,fdr_peptides))
  colnames(resulttable) <- c('SibPepCorr_cutoff' ,
                             'n_targetProteins' , 'n_true_targetProteins','n_decoyProteins' , 'proteinFDR',
                             'n_targetPeptides' , 'n_true_targetPeptides','n_decoyPeptides' , 'peptideFDR')
  }else if(fdr_type == "protein"){

    fdr_protein <- pmin(1, FFT*decoyProteins/(targetProteins))
    true_targetProteins <- as.integer(targetProteins * (1-fdr_protein))
    resulttable <- as.data.table(cbind(corrTest,
                                       targetProteins, true_targetProteins, decoyProteins,fdr_protein))
    colnames(resulttable) <- c('SibPepCorr_cutoff' ,
                               'n_targetProteins' , 'n_true_targetProteins','n_decoyProteins' , 'proteinFDR')

  }
  ## Only report sensible values
  resulttable <- resulttable[true_targetProteins >= 0.2*max(true_targetProteins)]

  # OUTPUT
  if (CSV == TRUE) {
    write.csv(resulttable,gsub("$|\\.csv$", ".csv", name), quote = FALSE, row.names = FALSE)
  }
  if (PDF == TRUE) {
    pdf(gsub("$|\\.pdf$", ".pdf", name))
  }

  if(plot){
    if(fdr_type =="protein"){
      xaxis <- resulttable$proteinFDR
      targets <- resulttable$n_targetProteins
      true_targets <- resulttable$n_true_targetProteins
    }else if(fdr_type =="peptide"){
      xaxis <- resulttable$peptideFDR
      targets <- resulttable$n_targetPeptides
      true_targets <- resulttable$n_true_targetPeptides
    }

    plot(xaxis, targets, type = "l", lty = 2,
       lwd = 2, main = paste("Precision Recall Curve"), xlim = c(0,0.1) ,xlab='FDR',ylab='n',
       ylim = c(0, 1.02* max(targets)))
    lines(xaxis, true_targets, lty = 1, lwd = 2)
    legend("bottomright", lty = c(2,1), legend = c("all target proteins", "true target proteins"))
  }
  if(PDF == TRUE){
    dev.off()
  }

  return(resulttable)
}
