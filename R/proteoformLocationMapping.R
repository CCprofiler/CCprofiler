#' Evaluate proteoform sequence proximity
#' @param traces Object of class traces or tracesList.
#' @return Traces with annotated sequence proximity stats.
#' @export
evaluateProteoformLocation <- function(traces, 
                                       adj.method = "fdr", 
                                       name="NormalizedSD", 
                                       minPepPerProtein = 6, 
                                       minPepPerProteoform = 3, ...){
  UseMethod("evaluateProteoformLocation", traces)
}

#' @describeIn evaluateProteoformLocation Evaluate proteoform sequence proximity
#' @export
evaluateProteoformLocation.traces <- function(traces, adj.method = "fdr", name="NormalizedSD", minPepPerProtein = 6, minPepPerProteoform = 3){
  traces$trace_annotation[, n_proteoforms := length(unique(proteoform_id)), by=c("protein_id")]
  traces$trace_annotation[, n_peptides := length(unique(id)), by=c("protein_id")]
  traces$trace_annotation[, n_peptides_per_proteoform := length(unique(id)), by=c("protein_id","proteoform_id")]
  medianPerProt <- traces$trace_annotation[, { 
    sub = unique(subset(.SD, select = c("proteoform_id","n_peptides_per_proteoform")))
    medianN = median(sub$n_peptides_per_proteoform)
    .(median_peptides_per_proteoform = as.double(medianN))}, by="protein_id"]
  # proteins_withoutProteoforms <- unique(traces$trace_annotation[n_proteoforms==1]$protein_id)
  proteins_withProteoforms <- unique(traces$trace_annotation[n_proteoforms!=1][n_peptides >= minPepPerProtein]$protein_id)
  proteins_withProteoforms <- proteins_withProteoforms[proteins_withProteoforms %in% medianPerProt[median_peptides_per_proteoform >= minPepPerProteoform]$protein_id]
  
  traces_toTest <- subset(traces, trace_subset_ids=proteins_withProteoforms, trace_subset_type="protein_id")
  testRandPep <- testRandomPeptides(traces_toTest)
  testRandPepStats <- plotRealVsRandom(testRandPep,score="NormalizedSD")
  
  testRandPepStats[, "genomLocation_pval_min" := min(genomLocation_pval), by="protein_id"]
  testRandPepStats_prot <- unique(subset(testRandPepStats, select=c("protein_id","genomLocation_pval_min")))
  
  testRandPepStats_prot$genomLocation_pval_adj_prot <- p.adjust(testRandPepStats_prot$genomLocation_pval_min, adj.method)
  
  pi0_prot <- pi0est(testRandPepStats_prot$genomLocation_pval_min, pi0.method = "bootstrap", na.rm=T)$pi0
  testRandPepStats_prot$pi0_prot <- pi0_prot
  
  qobj_prot <- qvalue(testRandPepStats_prot$genomLocation_pval_min, pi0=pi0_prot)
  testRandPepStats_prot$genomLocation_qval_prot <- qobj_prot$qvalues
  
  testRandPepStats <- merge(testRandPepStats, testRandPepStats_prot, by=c("protein_id","genomLocation_pval_min"))
  
  testRandPepStats$genomLocation_pval_adj <- p.adjust(testRandPepStats$genomLocation_pval, adj.method)
  
  pi0 <- pi0est(testRandPepStats$genomLocation_pval, pi0.method = "bootstrap", na.rm=T)$pi0
  qobj <- qvalue(testRandPepStats$genomLocation_pval, pi0=pi0)
  testRandPepStats$genomLocation_qval <- qobj$qvalues
  testRandPepStats$pi0 <- qobj$pi0
  
  pdf(paste0(name,".pdf"),width=3,height=3)
  p <- ggplot(testRandPepStats,aes(x=genomLocation_pval)) +
    geom_histogram(bins=25) +
    theme_classic()
  print(p)
  a <- ggplot(testRandPepStats,aes(x=genomLocation_pval_adj)) +
    geom_histogram(bins=25) +
    theme_classic()
  print(a)
  q <- ggplot(testRandPepStats,aes(x=genomLocation_qval)) +
    geom_histogram(bins=25) +
    theme_classic()
  print(q)
  p <- ggplot(testRandPepStats_prot,aes(x=genomLocation_pval_min)) +
    geom_histogram(bins=25) +
    theme_classic()
  print(p)
  a <- ggplot(testRandPepStats_prot,aes(x=genomLocation_pval_adj_prot)) +
    geom_histogram(bins=25) +
    theme_classic()
  print(a)
  q <- ggplot(testRandPepStats_prot,aes(x=genomLocation_qval_prot)) +
    geom_histogram(bins=25) +
    theme_classic()
  print(q)
  dev.off()
  traces$trace_annotation <- merge(traces$trace_annotation, testRandPepStats,
                                   all.x=T,all.y=F,by=c("protein_id","proteoform_id"),sort=F)
  .tracesTest(traces)
  return(traces)
}

#' @describeIn evaluateProteoformLocation Evaluate proteoform sequence proximity
#' @export
evaluateProteoformLocation.tracesList <- function(traces, 
                                                  adj.method = "fdr", 
                                                  name="NormalizedSD", 
                                                  minPepPerProtein = 6, 
                                                  minPepPerProteoform = 3){
  .tracesListTest(traces)
  
  traces_int <- integrateTraceIntensities(traces,
                                          design_matrix = NULL,
                                          integrate_within = NULL,
                                          aggr_fun = "sum")
  
  traces_int_ann <- evaluateProteoformLocation.traces(traces_int,
                                                      adj.method = adj.method, 
                                                      name=name, 
                                                      minPepPerProtein = minPepPerProtein, 
                                                      minPepPerProteoform = minPepPerProteoform)
  
  annTable = subset(traces_int_ann$trace_annotation, select=c("id","protein_id", 
                                                              names(traces_int_ann$trace_annotation)[!names(traces_int_ann$trace_annotation) %in% 
                                                                                                       names(traces[[1]]$trace_annotation)]))
  
  traces_ann <- lapply(traces, function(x){
    x$trace_annotation <- merge(x$trace_annotation,
                                annTable,by=c("protein_id","id"),sort=F,all.x=T,all.y=F)
    x
  })
  class(traces_ann) <- "tracesList"
  
  .tracesListTest(traces_ann)
  return(traces_ann)
}





testRandomPeptides <- function(traces){
  ann <- traces$trace_annotation
  proteins <- unique(ann$protein_id)
  res <- lapply(proteins,calculateScoresPerProtein,data=ann)
  names(res) <- proteins
  return(res)
}

calculateScoresPerProtein <- function(p,data,n_random=1000,seed=123){
  data_sub <- subset(data,protein_id==p)
  setkeyv(data_sub, c("protein_id","PeptidePositionStart"))
  data_sub[,PeptidePositionStartRank := seq_len(.N), by="protein_id"]
  n_peptides <- unique(data_sub$n_peptides)
  set.seed(seed)
  random_list <- lapply(c(1:n_random),function(i){sample(n_peptides, n_peptides, replace=F)})
  proteoforms <- unique(data_sub$proteoform_id)
  res <- lapply(proteoforms, calculateScoresPerProteoform, data=data_sub, random_list=random_list)
  names(res) <- proteoforms
  return(res)
}

calculateScoresPerProteoform <- function(p, data, random_list){
  idx <- which(data$proteoform_id == p)
  data_sub <- data[idx]
  rand_sub <- lapply(random_list,function(x){x[idx]})
  real_minDist <- getMinDist(data_sub$PeptidePositionStartRank)
  random_minDist <- lapply(rand_sub,getMinDist)
  real_n_notOne <- length(which(real_minDist!=1))
  random_n_notOne <- lapply(random_minDist,function(x){length(which(x!=1))})
  n_pepPerProteoform <- length(idx)
  real_f_notOne <- length(which(real_minDist!=1))/n_pepPerProteoform
  random_f_notOne <- lapply(random_minDist,function(x){length(which(x!=1))/n_pepPerProteoform})
  real_NormalizedSD <- getNormalizedSD(data_sub$PeptidePositionStartRank)
  random_NormalizedSD <- lapply(rand_sub,getNormalizedSD)
  list(
    n_notOne = list (
    real=real_n_notOne,
    random=random_n_notOne),
    f_notOne = list (
    real=real_f_notOne,
    random=random_f_notOne),
    NormalizedSD = list(
    real=real_NormalizedSD,
    random=random_NormalizedSD)
  )
}

getMinDist <- function(v){
  d <- dist(v,method="manhattan")
  m <- as.matrix(d)
  m[m == 0] <- Inf
  min_d <- apply(m, 1, min)
}

getNormalizedSD <- function(v){
  (sd(v)*sqrt((length(v)-1)/length(v)))/sqrt(((length(v)^2)-1)/12)}

plotRealVsRandom <- function(res,score="NormalizedSD"){
  proteins <- names(res)
  pdf(paste0(score,".pdf"),width=3,height=3)
  res <- lapply(proteins,plotRealVsRandomPerProt,res=res,score=score)
  dev.off()
  out <- do.call(rbind,res)
  return(out)
}

plotRealVsRandomPerProt <- function(protein,res,score){
  proteoforms <- names(lapply(res[[protein]],unlist))
  res <- lapply(proteoforms,plotRealVsRandomPerProteoform,protein=protein,res=res,score=score)
  out <- do.call(rbind,res)
  return(out)
}

plotRealVsRandomPerProteoform <- function(proteoform,protein,res,score){
  real <- res[[protein]][[proteoform]][[score]]$real
  random <- unlist(res[[protein]][[proteoform]][[score]]$random)
  #zscore <- (real-mean(random))/sd(random)
  #zcritical <- qnorm(0.05, mean = mean(random), sd = sd(random), lower.tail = TRUE, log.p = FALSE)
  n_rand_total <- length(random)
  n_rand_smallerReal <- length(which(random <= real))
  if (is.na(real)){
    p_rand <- 1
  } else {
    p_rand <- 1/n_rand_total*n_rand_smallerReal
  }
  dt <- data.table(random=random)
  if (!is.na(real)){
    p <- ggplot(dt,aes(x=random)) +
      geom_histogram(bins=30) +
      geom_vline(xintercept = real, colour="red") +
      theme_classic() +
      ggtitle(paste0(proteoform,
                     #"\n z-score = ",zscore,
                     #"\n z-critical = ",zcritical,
                     "\n p-rand = ",p_rand))
    print(p)
  }
  out <- data.table(protein_id=protein,proteoform_id=proteoform,genomLocation_pval=p_rand)
  return(out)
}

#' Plot peptide clusters by relative sequence position
#' @param traces Object of class traces or tracesList.
#' @param protein Character string of protein_id to plot.
#' @param PDF logical if PDF should written, default=FALSE.
#' @param closeGaps logical position gaps should be ignored, default=FALSE.
#' @return plot
#' @export
plotPeptideCluster <- function(traces,protein, PDF=FALSE, closeGaps=FALSE){
  traces <- subset(traces, protein, trace_subset_type="protein_id")
  if ("genomic_coord" %in% names(traces)) {
    traces$trace_annotation[,exon_id:=lapply(id,getGenomicCoord, traces=traces), by="id"]
    traces$trace_annotation[,min_exon_id_start:=min(PeptidePositionStart), by="exon_id"]
    getExonLevels <- unique(subset(traces$trace_annotation, select=c("exon_id","min_exon_id_start")))
    setorder(getExonLevels, min_exon_id_start)
    traces$trace_annotation$exon_id <- factor(traces$trace_annotation$exon_id, levels=getExonLevels$exon_id)
    getPalette = colorRampPalette(brewer.pal(9, "Spectral"))
  }
  dt <- subset(traces$trace_annotation,protein_id==protein)
  setkeyv(dt, c("protein_id","PeptidePositionStart"))
  dt[,PeptidePositionStartRank := seq_len(.N), by="protein_id"]
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#756bb1","#1c9099")
  dt$cluster <- as.factor(dt$cluster)
  if (PDF){
    if ("genomic_coord" %in% names(traces)) {
      pdf(paste0(protein,"_sequence_cluster.pdf"),width=10,height=5)
    } else {
      pdf(paste0(protein,"_sequence_cluster.pdf"),width=10,height=3)
    }
  }
  if (closeGaps) {
    q <- ggplot(dt,aes(x=PeptidePositionStartRank,
                       y=1,
                       fill=cluster)) +
      geom_bar(stat="identity") + theme_classic() +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title= element_blank(),
            axis.line = element_blank()) +
      theme(legend.position="bottom") +
      scale_fill_manual(values=cbPalette) 
    if ("genomic_coord" %in% names(traces)) {
      e <- ggplot(dt,aes(x=PeptidePositionStartRank,
                         y=1,
                         fill=exon_id)) +
        geom_bar(stat="identity") + theme_classic() +
        theme(axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title= element_blank(),
              axis.line = element_blank()) +
        theme(legend.position="top") +
        scale_fill_manual(values = getPalette(length(unique(dt$exon_id))))
      f <- ggarrange(e, q, 
                     labels = c("", ""),
                     ncol = 1, nrow = 2)
      print(annotate_figure(f, fig.lab = paste0(protein," : ",unique(dt$Gene_names))))
    } else {
      print(q + ggtitle(paste0(protein," : ",unique(dt$Gene_names)))) 
    }
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
      scale_fill_manual(values=cbPalette) 
    if ("genomic_coord" %in% names(traces)) {
      e <- ggplot(dt,aes(x=PeptidePositionStart,
                         y=1,
                         fill=exon_id)) +
        geom_bar(stat="identity") + theme_classic() +
        theme(axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title= element_blank(),
              axis.line = element_blank()) +
        theme(legend.position="top") +
        scale_fill_manual(values = getPalette(length(unique(dt$exon_id))))
      f <- ggarrange(e, q, 
                     labels = c("", ""),
                     ncol = 1, nrow = 2)
      print(annotate_figure(f, fig.lab = paste0(protein," : ",unique(dt$Gene_names))))
    } else {
      print(q + ggtitle(paste0(protein," : ",unique(dt$Gene_names)))) 
    }
  }
  if (PDF){
    dev.off()
  }
}
