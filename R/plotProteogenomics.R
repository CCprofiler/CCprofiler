#' Plot peptides in their genomic context
#' @param traces Object of class traces with an additional item =genomic_coord=
#' @param gene_id character, the gene_id of the protein to be plotted
#' @param ensdb EnsDb object or path to EnsDb database. If not NULL, all transcripts
#' for the specified gene_id will be plotted.
#' @param colorMap named character vector containing valid color specifications for plotting.
#' The names of the vector must correspond to the ids of the peptides to be plotted.
#' @param plot boolean, whether to plot or to return the plot.
#' @import ggbio
#' @export

plotPeptidesInGenome <- function(traces,
                                 gene_id,
                                 ensdb=NULL,
                                 colorMap=NULL,
                                 highlight=NULL,
                                 plot=T){

  if(!is.null(ensdb)){
    edb <- validateEnsDB(ensdb)
    p_gm <- ggbio::autoplot(edb, GeneIdFilter(gene_id))
  }

  tr_sub <- subset.traces(traces, gene_id, "gene_id")

  if(is.null(tr_sub[["genomic_coord"]])){
    message("No genomic Coordinates found. Trying to fetch them...")
    tr_sub <- annotateGenomicCoordinates(tr_sub, db=edb)
  }
  pepRanges <- tr_sub$genomic_coord

  gr <- unlist(pepRanges)
  mcols(gr)$id <- names(gr)
  pepRanges <- relist(gr, pepRanges)


  ## Annotate which traces to highlight
  if(!is.null(highlight)){
    mcols(pepRanges)$outlier <- gsub("\\(.*?\\)","",names(pepRanges)) %in% gsub("\\(.*?\\)","",highlight)
    mcols(pepRanges)$outlier <- sapply(mcols(pepRanges)$outlier, ifelse, "a","b")
    if(!any(mcols(pepRanges)$outlier == "a")) highlight <- NULL
  }else{
    mcols(pepRanges)$outlier <- "a"
  }

  ## Create a reproducible coloring for the peptides plotted
  if(!is.null(colorMap)){
    if(!all(unique(names(pepRanges)) %in% names(colorMap))){
      stop("Invalid colorMap specified. Not all traces to be plotted are contained in the colorMap")
    }
  }else{
    colorMap <- createGGplotColMap(names(pepRanges))
  }

  p_pep <- ggplot() +
    geom_alignment(pepRanges, aes(fill=id, alpha=outlier),
                   type="model", group.selfish=F) +
    theme(legend.position="none") +
    scale_fill_manual(values=colorMap) +
    scale_alpha_manual(values = c(a=1, b=0.3))

  if(!is.null(ensdb)){
    tr <- tracks(GeneModel=p_gm, Peptides=p_pep, heights = c(2,1))
  }else{
    tr <- tracks(Peptides=p_pep)
  }
  p_tr <- plot(tr_sub, plot=F, highlight=highlight, colorMap=colorMap)
  tr <- as(tr, "grob")
  if(plot){
    grid.arrange(p_tr, tr, ncol=1)
  }else{
    return(grid.arrange(p_tr, tr, ncol=1))
  }
}
