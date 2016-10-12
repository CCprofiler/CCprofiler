ROC.SibPepCorr <- function(Traces, FFT = 0.50, Stepsize =0.001 , Summary = FALSE, PDF = FALSE) {
  
  trace_annotation <- Traces$trace_annotation
  trace_annotation$DECOY <- 0
  trace_annotation$DECOY[grep('DECOY', trace_annotation$protein_id)] <- 1
  corr_range <- sort(unique(trace_annotation$SibPepCorr))
  corr.test <- seq(min(corr_range), max(corr_range), Stepsize)
  ncorr.test <- length(corr.test)
  target_proteins <- numeric(length = ncorr.test)
  decoy_proteins <- numeric(length = ncorr.test)
  target_protein_ids <- character(length = ncorr.test)
  
  for (i in 1:ncorr.test) {
    
    targetprots <- unique(trace_annotation[SibPepCorr >= corr.test[i] & DECOY == 0]$protein_id)
    target_protein_ids[i] <- paste(targetprots, collapse = ",")
    target_proteins[i] <- length(targetprots)
    decoyprots <- unique(trace_annotation[SibPepCorr >= corr.test[i] & DECOY == 1]$protein_id)
    decoy_proteins[i] <- length(decoyprots)
    
    }
  
  fdr_protein <- FFT*decoy_proteins/(target_proteins+decoy_proteins)
  true_target_proteins <- as.integer(target_proteins * (1-fdr_protein))
  resulttable <- as.data.table(cbind(corr.test, target_proteins, true_target_proteins, decoy_proteins,fdr_protein))
  colnames(resulttable) <- c('Tested Corr' , 'Target proteins' , 'True target proteins','Decoys' , 'FDR')
  
  if (Summary == TRUE) {
    write.table(resulttable,'FDR_based_on_SibPepCorr.txt',quote = FALSE,sep = '\t',row.names = FALSE)
  }
  
  if(PDF == TRUE){
    pdf("CorrFilter_SibPepCorr_ROC.pdf")
  }
  plot(resulttable$FDR, resulttable$`Target proteins`, type = "l", lty = 2, main = paste("ROC Curve"), xlim = c(0,0.1) ,xlab='FDR',ylab='Proteins')
  lines(resulttable$FDR, resulttable$`True target proteins`, lty = 1)
  legend("bottomright", lty = c(2,1), legend = c("all target proteins", "true target proteins"))
  if(PDF == TRUE){
    dev.off()
  }
}
  
  
  
  
  
  
  
  
  