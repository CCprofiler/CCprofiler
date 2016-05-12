#' From protein level Traces, generate a Network table with pairwise profile
#' correlation score
#' 
#' @param Traces A Traces object containing wide format elution profiles.
#' @param correlation_cutoff Minimum pearson correlation value to report an
#'     edge between two proteins.
#' @param write_file Logical (TRUE/FALSE) whether an output file named
#'     [systime]_NetworkTable.csv shall be written to the working folder.
#' @return A data.table containing columns UniProt_A UniProt_B and
#'     pearson_correlation, readable e.g. by Cytoscape.
#' @export
generateNetworkTable <- function(Traces, correlation_cutoff=0.9, write_csv=TRUE) {
    data <- Traces$protein.traces
    # calculate protein level correlation matrix
    data.traces <- as.matrix(data[,2:ncol(data)])
    rownames(data.traces) <- data$protein_id
    corrmtrx <- cor(t(data.traces))
    corrlist <- corrmtrx
    corrlist[lower.tri(corrlist,diag=TRUE)] <- NA
    corrlist <- as.data.frame(as.table(corrlist))
    corrlist <- na.omit(corrlist)
    corrlist <- corrlist[order(-abs(corrlist$Freq)), ]
    names(corrlist) <- c('Uniprot_A', 'Uniprot_B', 'pearson_correlation')
    corrlist.s <- subset(corrlist, pearson_correlation >= correlation_cutoff)
    if (write_csv == TRUE) {
        write.csv(corrlist.s, file=paste0(deparse(substitute(Traces)),
                                          '_NetworkTable.csv', quote=FALSE,
                                          row.names=FALSE))
    }
    return(corrlist.s)
}
