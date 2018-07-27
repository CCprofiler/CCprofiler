#' Plot peptides in their genomic context
#' @param traces Object of class traces with an additional item =genomic_coord=
#' @param geneId character, the geneId of the protein to be plotted
#' @param ensdb EnsDb object or path to EnsDb database. If not NULL, all transcripts
#' for the specified geneId will be plotted.
#' @param colorMap named character vector containing valid color specifications for plotting.
#' The names of the vector must correspond to the ids of the peptides to be plotted.
#' @param plot boolean, whether to plot or to return the plot.
#' @param ... Additional arguments to be passed to plot.traces
#' @import ggbio
#' @export

plotPeptidesInGenome <- function(traces,
                                 geneId,
                                 ensdb=NULL,
                                 colorMap=NULL,
                                 highlight=NULL,
                                 highligt_detected_transcripts=TRUE,
                                 ...){


  tr_sub <- subset.traces(traces, geneId, "gene_id")

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
    theme_bw() +
    theme(legend.position="none") +
    scale_fill_manual(values=colorMap) +
    scale_alpha_manual(values = c(a=1, b=0.3))

  tr <- tracks(Peptides=p_pep, heights=1)

  if(!is.null(highlight)){
    p_pepO <- ggplot() + geom_alignment(pepRanges[mcols(pepRanges)$outlier == "a",], aes(fill = id),
                                        color = NA, type = "model", group.selfish = F) + 
      theme_bw() +
      theme(legend.position = "none") +
      scale_fill_manual(values = colorMap)
    h <- tr@heights
    tr <- tr + tracks(OutlPeptides = p_pepO)
    tr@heights <- c(h,1)
  }

  ## Plot the gene model
  if(!is.null(ensdb)){
    edb <- validateEnsDB(ensdb)
    ## get transcripts to highlight
    hl_transcripts <- NULL
    if(highligt_detected_transcripts){
      ensProtIds <- getEnsemblProteins(tr_sub)
      ensTxIds <- ensembldb::select(ensdb, keys = ~protein_id == ensProtIds, columns = c("TXID"))$TXID
      hl_transcripts <- ~ tx_id == ensTxIds
    }
    p_gm <- plotGeneModel(edb, geneId, highlight = hl_transcripts)
    if(class(p_gm) != "try-error"){
      h <- tr@heights
      tr <- tracks(Gene_model = p_gm) + tr
      tr@heights <- c(4, h)
    }
  }
  p_tr <- plot(tr_sub, plot=F, highlight=highlight, colorMap=colorMap, ...)
  tr <- as(tr, "grob")
  return(grid.arrange(p_tr, tr, ncol=1))
}


#' Plot the known transcripts for a gene of interest
#' @param ensdb EnsDb object or path to EnsDb database. If not NULL, all transcripts
#' for the specified geneId will be plotted.
#' @param geneId character, the geneId of the protein to be plotted
#' @param filter formula, Filtering expression to be passed to AnnotationFilter
#' @importFrom ensembldb select
#' @export

plotGeneModel <- function(ensdb, geneId,
                          filter=~ gene_id == geneId & tx_biotype == "protein_coding",
                          highlight=NULL){
  ## TODO explore truncate.gaps to make the exons bigger
  ## see https://bioconductor.org/packages/release/bioc/vignettes/biovizBase/inst/doc/intro.pdf
  fillMap <- c(cds="white", exon="black", utr="blue")
  alphaMap <- c(a=1, b=0.5)

  p_gm <- try({
    gene_filt <- AnnotationFilter(filter)
    p <- ggbio::autoplot(ensdb, gene_filt, aes(fill=type),
                         cds.rect.h=0.25,
                         exon.rect.h=0.25,
                         utr.rect.h=0.25) +
      theme_bw() +
      theme(legend.position="none") +
      scale_fill_manual(values = fillMap)

    if(!is.null(highlight)){
      hl_ids <- ensembldb::select(ensdb, keys = highlight, columns = c("ENTREZID","TXID"))$TXID

      ## ggplot uses the .GlobalEnv for aes evaluation by default, here we want to use a local Env
      ## Pass all external variables via eval(substitute)
      eval(substitute(
        expr =
          p <- p + aes(alpha=ifelse(tx_id %in% hl_ids, "a", "b")) +
            scale_alpha_manual(values = alphaMap),
        env = list(hl_ids=hl_ids)))
    }
    p
  })
  return(p_gm)
}
