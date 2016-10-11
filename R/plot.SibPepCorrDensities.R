
# For some pretty weird reason
.datatable.aware=TRUE




plot.SibPepCorrDensities <- function(Traces, PDF = FALSE){
  
  if(PDF == TRUE){
    pdf("CorrFilter_SibPepCorr_distributions_target_decoy.pdf")
  }
  
  
  trace_annot <- Traces$trace_annotation
  dens_all <- density(na.omit(trace_annot$SibPepCorr))
  dens_targets <- density(na.omit(trace_annot[grep("^1/", trace_annot$ProteinName),]$SibPepCorr))
  dens_decoys <- density(na.omit(trace_annot[grep("^DECOY_1/", trace_annot$ProteinName),]$SibPepCorr))
  
  plot(dens_targets$x, dens_targets$y*dens_targets$n, lty = 1, lwd = 3,
       , type = "l", ylab = "scaled frequency", xlab = "SibPepCorr",
       main = "SEC-SWATH\n Sibling Peptide Correlation target vs decoy peptides")
  
  lines(dens_decoys$x, dens_decoys$y*dens_decoys$n, lty = 2, lwd = 3)
  
  legend("topleft", lty = c(1,3), lwd = c(3,3), legend = c("target peptides", "decoy peptides"))
  
  
  if(PDF == TRUE){
    dev.off()
  }
}

