#' getGlobalPairwiseCorr
#' @description From peptide or protein level Traces, generate a Network table with pairwise profile
#' correlation score
#' @param traces An object of type \code{traces}.
#' @param correlation_cutoff Minimum pearson correlation value to report an edge between two proteins, default=0.9.
#' @param CSV Logical (TRUE/FALSE) whether an output file named GlobalPairwiseCorrelationTable.csv shall be written to the working folder.
#' @return A data.table containing columns peptide/protein_A peptide/protein_B and pearson_correlation, readable e.g. by Cytoscape.
#' @export
getGlobalPairwiseCorr <- function(traces, correlation_cutoff = 0.9, write_csv = TRUE) {
	data.traces=getIntensityMatrix(traces)
	corrmtrx<-cor(t(data.traces))
	corrlist<-corrmtrx
	corrlist[lower.tri(corrlist,diag=TRUE)]=NA
	corrlist=as.data.frame(as.table(corrlist))
	corrlist=na.omit(corrlist)
	corrlist=corrlist[order(-abs(corrlist$Freq)),]
	if(traces$trace_type=="peptide") {
		names(corrlist)<-c("peptide_A", "peptide_B", "pearson_correlation")
	} else if (traces$trace_type=="protein") {
		names(corrlist)<-c("protein_A", "protein_B", "pearson_correlation")
	} else {
		stop("Invalid trace type. Traces should be of trace_type peptide or protein.")
	}
	corrlist.s <- subset(corrlist, pearson_correlation >= correlation_cutoff)
	if (write_csv == TRUE) {
		write.csv(corrlist.s, file="GlobalPairwiseCorrelationTable.csv", quote = FALSE, row.names = FALSE)
	}
	return(corrlist.s)
}
