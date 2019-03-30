
getGenomicCoord <- function(id,traces){
  x <- traces$genomic_coord[[id]]
  return(x$exon_id[1])
}

countMinSwaps <- function(dt){
  if (length(unique(dt$proteoform_id)) == 1) {
    return(0)
  } else {
    max_table <- max(table(dt$proteoform_id))
    minimum_swaps <- nrow(dt)-max_table
    return(minimum_swaps)
  }
}

evaluateExonSupport <- function(traces,n_random=1000,seed=123){
  traces$trace_annotation[,exon_id:=lapply(id,getGenomicCoord, traces=traces), by="id"]
  traces$trace_annotation[,n_exons:=length(unique(exon_id)), by="protein_id"]
  traces$trace_annotation[,n_peptides_per_exon := .N, by=c("exon_id","protein_id")]
  traces$trace_annotation[,min_peptides_per_exon := min(n_peptides_per_exon), by="protein_id"]
  traces$trace_annotation[,max_peptides_per_exon := max(n_peptides_per_exon), by="protein_id"]
  medianPerProt <- traces$trace_annotation[, { 
    sub = unique(subset(.SD, select = c("exon_id","n_peptides_per_exon")))
    medianN = median(sub$n_peptides_per_exon)
    .(median_peptides_per_exon = as.double(medianN))}, by="protein_id"]
  #proteins <- unique(traces$trace_annotation[max_peptides_per_exon>2]$protein_id)
  # only use proteins with a medium >= 2 peptides per exon
  proteins <- unique(medianPerProt[median_peptides_per_exon >= 2]$protein_id)
  # only use genes with > 1 exon
  proteins <- proteins[proteins %in% traces$trace_annotation[n_exons > 1]$protein_id]
  res_prot <- lapply(proteins, function(p){
    traces_sub <- subset(traces,trace_subset_ids=p,trace_subset_type="protein_id")
    set.seed(seed)
    all_exons <- traces_sub$trace_annotation$exon_id
    random_exon_list <- lapply(c(1:n_random),function(i){sample(all_exons, length(all_exons), replace=F)})
    random_exon_list <- lapply(random_exon_list,function(l){
      data.table(exon_id=l,proteoform_id=traces_sub$trace_annotation$proteoform_id)
    })
    exons <- unique(all_exons)
    res_exon <- lapply(exons, function(e){
      traces_exon <- subset(traces,trace_subset_ids=e,trace_subset_type="exon_id")
      traces_exon_dt <- subset(traces_exon$trace_annotation,select=c("exon_id","proteoform_id"))
      random_exon_list_sub <- lapply(random_exon_list,function(l){
        subset(l,exon_id==e)
      })
      exon_observed <- countMinSwaps(traces_exon_dt)
      exon_random <- unlist(lapply(random_exon_list_sub,countMinSwaps))
      return(list(observed=exon_observed,random=exon_random))
    })
    names(res_exon) <- exons
    return(res_exon)
  })
  names(res_prot) <- proteins
  final_res <- lapply(res_prot,function(l){
    sum_observed <- sum(unlist(lapply(l, `[[`, 1)))
    sum_random <- rowSums(as.data.table(lapply(l, `[[`, 2)))
    nExons <- length(l)
    return(list(observed=sum_observed,random=sum_random,nExons=nExons))
  })
  return(final_res)
}

plotRealVsRandomSwaps <- function(res){
  proteins <- names(res)
  pdf(paste0("RandomSwaps.pdf"),width=3,height=3)
  res <- lapply(proteins,plotRealVsRandomPerProtein,res=res)
  dev.off()
  out <- do.call(rbind,res)
  return(out)
}

plotRealVsRandomPerProtein <- function(protein,res){
  real <- res[[protein]]$observed
  random <- res[[protein]]$random
  n_rand_total <- length(random)
  n_rand_smallerReal <- length(which(random <= real))
  n_min_mistake <- length(which(random == 0))
  p_rand <- 1/n_rand_total*n_rand_smallerReal
  min_possible_pval <- 1/n_rand_total*n_min_mistake
  nExons <- res[[protein]]$nExons
  dt <- data.table(random=random)
  p <- ggplot(dt,aes(x=random)) +
    geom_histogram(binwidth = 1) +
    geom_vline(xintercept = real, colour="red") +
    theme_classic() +
    ggtitle(paste0(protein,
      "\n p-rand = ",p_rand))
  print(p)
  out <- data.table(protein_id=protein,exon_pval=p_rand,nExonMistakes=real,nExons=nExons,min_possible_pval = min_possible_pval)
  return(out)
}

#' Evaluate if proteoforms agree with genomic information of exons
#' @param traces Object of class traces or tracesList.
#' @param adj.method Character string which p-value adjustment to use.
#' Default "fdr".
#' @return Traces with exon evidence p-value per gene
#' @export
evaluateExonLocation <- function(traces, adj.method = "fdr", optional_filter = F){
  traces$trace_annotation[, n_proteoforms := length(unique(proteoform_id)), by=c("protein_id")]
  proteins_withoutProteoforms <- unique(traces$trace_annotation[n_proteoforms==1]$protein_id)
  proteins_withProteoforms <- unique(traces$trace_annotation[n_proteoforms!=1]$protein_id)
  traces_toTest <- subset(traces, trace_subset_ids=proteins_withProteoforms, trace_subset_type="protein_id")

  exonRes <- evaluateExonSupport(traces_toTest)
  exonStats <- plotRealVsRandomSwaps(exonRes)
  # remove results for proteins with nExons < n_peptides-1
  nPeps <- unique(subset(traces$trace_annotation,select=c("protein_id","n_peptides")))
  exonStats <- merge(exonStats, nPeps,
    all.x=T,all.y=F,by=c("protein_id"),sort=F)
  #exonStats <- subset(exonStats, nExons < n_peptides-1)
  exonStats$exon_pval_adj <- p.adjust(exonStats$exon_pval, adj.method)
  exonStats$min_possible_pval_adj <- p.adjust(exonStats$min_possible_pval, adj.method)
  exonStats[,n_peptides:=NULL]
  # filter proteins with exon_pvals == 1 or (optional) remove proteins
  # with low powers as well as exon_pvals == 1
  if (optional_filter == FALSE) {
    exonStats <- exonStats
    # exonStats <- subset(exonStats, exon_pval != 1)
  } else {
    exonStats <- subset(exonStats, min_possible_pval < 0.05 | exon_pval != 1)
  }

  pdf("RealVsRandomSwaps_hist.pdf",width=3,height=3)
    p <- ggplot(exonStats,aes(x=exon_pval)) +
    geom_histogram(bins=25)+
    theme_classic()
    print(p)
    q <- ggplot(exonStats,aes(x=exon_pval_adj)) +
    geom_histogram(bins=25)+
    theme_classic()
    print(q)
  dev.off()

  pdf("min_possible_pval_hist.pdf",width=3,height=3)
    r <- ggplot(exonStats,aes(x=min_possible_pval)) +
    geom_histogram(bins=25)+
    theme_classic()
    print(r)
    s <- ggplot(exonStats,aes(x=min_possible_pval_adj)) +
    geom_histogram(bins=25)+
    theme_classic()
    print(s)
  dev.off()

  traces$trace_annotation <- merge(traces$trace_annotation, exonStats,
    all.x=T,all.y=F,by=c("protein_id"),sort=F)
  .tracesTest(traces)
  return(traces)
}

#' Plot peptide clusters by relative sequence position
#' @param traces Object of class traces or tracesList.
#' @param protein Character string of protein_id to plot.
#' @param PDF logical if PDF should written, default=FALSE.
#' @param closeGaps logical position gaps should be ignored, default=FALSE.
#' @return plot
#' @export
plotPeptideCluster <- function(traces,protein, PDF=FALSE, closeGaps=FALSE){
  dt <- subset(traces$trace_annotation,protein_id==protein)
  setkeyv(dt, c("protein_id","PeptidePositionStart"))
  dt[,PeptidePositionStartRank := seq_len(.N), by="protein_id"]
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#756bb1","#1c9099")
  dt$cluster <- as.factor(dt$cluster)
  if (PDF){
    pdf(paste0(protein,"_sequence_cluster.pdf"),width=10,height=3)
  }
  if (closeGaps) {
    p <- ggplot(dt,aes(x=PeptidePositionStartRank,
      y=1,
      fill=cluster)) +
      geom_bar(stat="identity") + theme_classic() +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title= element_blank(),
            axis.line = element_blank()) +
      theme(legend.position="bottom") +
      scale_fill_manual(values=cbPalette) +
      ggtitle(paste0(protein," : ",unique(dt$Gene_names)))
    print(p)
  } else {
    q <- ggplot(dt,aes(x=PeptidePositionStart,
      y=1,
      fill=cluster)) +
      geom_bar(stat="identity") + theme_classic() +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title= element_blank(),
            axis.line = element_blank()) +
      theme(legend.position="bottom") +
      scale_fill_manual(values=cbPalette) +
      ggtitle(paste0(protein," : ",unique(dt$Gene_names)))
    print(q)
  }
  if (PDF){
    dev.off()
  }
}
