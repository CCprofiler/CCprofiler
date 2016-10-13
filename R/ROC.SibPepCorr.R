ROC.SibPepCorr <- function(Traces,
                           FFT = 1,
                           Stepsize =0.001,
                           plot = FALSE,
                           PDF = FALSE,
                           CSV = FALSE) {
  
  trace_annotation <- Traces$trace_annotation
  trace_annotation$DECOY <- 0
  trace_annotation$DECOY[grep('^DECOY_', trace_annotation$protein_id)] <- 1
  corr_range <- sort(unique(trace_annotation$SibPepCorr))
  corr.test <- seq(min(corr_range), max(corr_range), Stepsize)
  ncorr.test <- length(corr.test)
  target_proteins <- numeric(length = ncorr.test)
  decoy_proteins <- numeric(length = ncorr.test)
  target_protein_ids <- character(length = ncorr.test)
  
  for (i in 1:ncorr.test) {
    
    targetprots <- unique(trace_annotation[SibPepCorr >= corr.test[i] & DECOY == 0]$protein_id)
    # target_protein_ids[i] <- paste(targetprots, collapse = ",")
    target_proteins[i] <- length(targetprots)
    decoyprots <- unique(trace_annotation[SibPepCorr >= corr.test[i] & DECOY == 1]$protein_id)
    decoy_proteins[i] <- length(decoyprots)
    
  }
  
  fdr_protein <- FFT*decoy_proteins/(target_proteins+decoy_proteins)
  true_target_proteins <- as.integer(target_proteins * (1-fdr_protein))
  resulttable <- as.data.table(cbind(corr.test, target_proteins, true_target_proteins, decoy_proteins,fdr_protein))
  colnames(resulttable) <- c('SibPepCorr_cutoff' , 'n_target_proteins' , 'n_true_target_proteins','n_decoys' , 'FDR')
  resulttable <- resulttable[true_target_proteins >= 0.2*max(true_target_proteins)]
  
  # OUTPUT
  if (CSV == TRUE) {
    write.csv(resulttable,'SibPepCorr_ROC.csv', quote = FALSE, row.names = FALSE)
  }
  if (PDF == TRUE) {
    pdf("SibPepCorr_ROC.pdf")
  }
  
  if(plot){
    plot(resulttable$FDR, resulttable$n_target_proteins, type = "l", lty = 2,
       lwd = 2, main = paste("ROC Curve"), xlim = c(0,0.1) ,xlab='FDR',ylab='n',
       ylim = c(0, 1.02* max(resulttable$n_target_proteins)))
    lines(resulttable$FDR, resulttable$n_true_target_proteins, lty = 1, lwd = 2)
    legend("bottomright", lty = c(2,1), legend = c("all target proteins", "true target proteins"))
  }
  if(PDF == TRUE){
    dev.off()
  }
  
  return(resulttable)
}