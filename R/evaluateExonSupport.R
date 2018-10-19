
getGenomicCoord <- function(id){
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
  traces$trace_annotation[,exon_id:=lapply(id,getGenomicCoord), by="id"]
  proteins <- unique(traces$trace_annotation$protein_id)
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
    return(list(observed=sum_observed,random=sum_random))
  })
  return(final_res)
}

plotRealVsRandomSwaps <- function(res){
  proteins <- names(res)
  pdf(paste0("RandomSwaps.pdf"))
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
  p_rand <- 1/n_rand_total*n_rand_smallerReal
  dt <- data.table(random=random)
  p <- ggplot(dt,aes(x=random)) +
    geom_histogram(bins=30) +
    geom_vline(xintercept = real, colour="red") +
    theme_classic() +
    ggtitle(paste0(protein,
      "\n p-rand = ",p_rand))
  print(p)
  out <- data.table(protein_id=protein,exon_pval=p_rand)
  return(out)
}
