#' Filter peptides for having
#' at least one high correlating sibling peptide
#' @param traces Object of class traces.
#' @param cutoff Numeric between 0 and 1. Minimum correlation of a peptide
#' width any sibling peptide.
#' @param plot logical,wether to print SibPepCorr density plot to R console.
#' Deafult is \code{FALSE}.
#' @param PDF logical, wether to print SibPepCorr density plot to a PDF file.
#' Deafult is \code{FALSE}.
#' @param name Character string with name of the plot, only used if
#' '\code{PDF=TRUE}.PDF file is saved under name.pdf. Default is "maxCorrHist".
#' @return Gene peptide list
#' @export
getGenePepList <- function(traces){
  genes <- unique(traces$trace_annotation$protein_id)
  intMat <- getIntensityMatrix(traces)
  genePepList <- lapply(genes, FUN = function(gene){
    peps <- traces$trace_annotation[protein_id == gene, id]
    res <- intMat[peps,]
    return(t(res))
  })
  names(genePepList) <- genes
  return(genePepList)
}

#' Filter peptides for having
#' at least one high correlating sibling peptide
#' @param traces Object of class traces.
#' @param cutoff Numeric between 0 and 1. Minimum correlation of a peptide
#' width any sibling peptide.
#' @param plot logical,wether to print SibPepCorr density plot to R console.
#' Deafult is \code{FALSE}.
#' @param PDF logical, wether to print SibPepCorr density plot to a PDF file.
#' Deafult is \code{FALSE}.
#' @param name Character string with name of the plot, only used if
#' '\code{PDF=TRUE}.PDF file is saved under name.pdf. Default is "maxCorrHist".
#' @param acrossConditions logical, if traces across replicates and conditions
#' should be integrated for filtering to ensure minimal data loss.
#' Deafult is \code{TRUE}.
#' @return Object of class traces filtered for peptide correlation.
#' @importFrom matrixStats rowMaxs
#' @export

filterByMaxCorr <- function(traces, cutoff = 0.85,
                            plot = FALSE, PDF=FALSE, name="maxCorrHist", ...) {
  genePeptideList <- getGenePepList(traces)
  maxCorrMatrices <- lapply(genePeptideList, function(gene){
    genecorr <- cor(gene)
    genecorr[genecorr == 1] <- NA
    maxcorr <- matrixStats::rowMaxs(genecorr, na.rm = T)
    names(maxcorr) <- rownames(genecorr)
    return(maxcorr)
  })

  names(maxCorrMatrices) <- names(genePeptideList)
  allMaxCorrs <- unlist(maxCorrMatrices)
  ## When unlisting protein and peptide names are concatenated with a .
  ## We need to be careful with peptides starting with . (e.g. .(UniMod))
  names(allMaxCorrs) <- gsub("^.*?\\.", "", names(allMaxCorrs))
  if (plot == TRUE){
    if (PDF) {
      pdf(paste0(name,".pdf"),width=3,height=3)
    }
    p <- ggplot(data.table(maxCorr=allMaxCorrs),
                aes(x=maxCorr)) +
                geom_histogram(bins=30) +
                theme_classic()
    plot(p)
    if (PDF) {
      dev.off()
    }
  }

  filterpeps <- names(allMaxCorrs)[allMaxCorrs > cutoff]

  traces_filt <- subset(traces, filterpeps, "id")
  return(traces_filt)
}

#' Combine tracesList across conditions to one traces object with all
#' the traces of all conditions added as additional fractions for
#' peptide clustering across conditions.
#' @param tracesList Object of class tracesList.
#' @return Object of class traces
#' @export
combineTracesMutiCond <- function(tracesList){
  .tracesListTest(tracesList)
  cond <- names(tracesList)
  idx <- seq_along(cond)
  traces <- lapply(idx, function(i){
    t <- tracesList[[i]]
    trac <- t$traces
    trac.m <- melt(trac,id.vars="id")
    trac.m[, cond := cond[[i]]]
    return(trac.m)
    })
  traces_combi <- rbindlist(traces)
  traces_all <- dcast(traces_combi,id ~ cond+variable,value.var="value",fill=0)
  names(traces_all) <- c("id",seq(1,ncol(traces_all)-1,1))
  setcolorder(traces_all, c(seq(1,ncol(traces_all)-1,1),"id"))

  trace_type_all <- tracesList[[1]]$trace_type

  trace_annotation <- lapply(tracesList, function(t){
    # Use the first element as a template annotation
    cols <- names(tracesList[[1]]$trace_annotation)
    cols <- cols[which(!cols %in% c("SibPepCorr","RepPepCorr","n_peptides"))]
    res <- subset(t$trace_annotation,select=cols)
    return(res)
    })
  trace_annotation_combi <- rbindlist(trace_annotation)
  trace_annotation_combi <- unique(trace_annotation_combi)
  trace_annotation_combi <- trace_annotation_combi[order(match(id,
                                                          traces_all$id))]

  fraction_annotation_all <- data.table(id=seq(1,ncol(traces_all)-1,1))

  combi_traces <- list(
    traces = traces_all,
    trace_type = trace_type_all,
    trace_annotation = trace_annotation_combi,
    fraction_annotation = fraction_annotation_all
  )

  if(!is.null(tracesList[[1]][["genomic_coord"]])){
	coo <- tracesList[[1]][["genomic_coord"]]
	for (i in c(2:length(tracesList))) {
		coo <- append(coo, tracesList[[i]][["genomic_coord"]], after = length(coo))
	}
	coo <- coo[!duplicated(names(coo))]
	combi_traces$genomic_coord <- coo
  }

  class(combi_traces) <- "traces"
  .tracesTest(combi_traces)
  return(combi_traces)
}
