
# For some pretty weird reason
.datatable.aware=TRUE


ROC.SibPepCorr <- function(Traces, FFT = 0.50, Stepsize =0.001 , Summary = FALSE) {
  
  trace_annot <- Traces$traces_annotation
  trace_annot$DECOY <- 0
  handle <- (substr(trace_annot$ProteinName,1,5) == 'DECOY')
  trace_annot$DECOY[handle] <- 1
  corr_range <- sort(unique(trace_annot$SibPepCorr))
  corr.test <- seq(min(corr_range), max(corr_range), Stepsize)
  
  ncorr.test <- length(corr.test)
  target_proteins <- numeric(length = ncorr.test)
  decoy_proteins <- numeric(length = ncorr.test)
  target_protein_ids <- character(length = ncorr.test)
  
  for (i in 1:ncorr.test) {
    
    targetprots <- unique(trace_annot[trace_annot$SibPepCorr >= corr.test[i] & trace_annot$DECOY == 0 ,]$ProteinName)
    target_protein_ids[i] <- paste(targetprots, collapse = ",")
    target_proteins[i] <- length(targetprots)
    decoyprots <- unique(trace_annot[trace_annot$SibPepCorr >= corr.test[i] & trace_annot$DECOY == 1 ,]$ProteinName)
    decoy_proteins[i] <- length(decoyprots)
    
    }
  
  fdr_protein <- FFT*decoy_proteins/(target_proteins+decoy_proteins)
  true_target_proteins <- as.integer(target_proteins * (1-fdr_protein))
  
  resulttable <- as.data.frame(cbind(corr.test, target_proteins, true_target_proteins, decoy_proteins,fdr_protein))
  colnames(resulttable) <- c('Tested Corr' , 'Target proteins' , 'True target proteins','Decoys' , 'FDR')
  
  if (Summary == TRUE) {
    write.table(resulttable,'FDR_based_on_SibPepCorr.txt',quote = FALSE,sep = '\t',row.names = FALSE)
  }
  
  
  plot(resulttable$FDR, resulttable$`Target proteins`, type = "l", lty = 2, main = paste("ROC Curve"), xlim = c(0,0.1) ,xlab='FDR',ylab='Proteins')
  
  lines(resulttable$FDR, resulttable$`True target proteins`, lty = 1)
  
  legend("bottomright", lty = c(2,1), legend = c("all target proteins", "true target proteins"))
  
}
  
  
  
  
  
  
  
  
  