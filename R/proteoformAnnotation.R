#' Extract matrix of trace intensities for every protein id
#' @param traces Object of class traces.
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

#' Calculate pairwise correlations of peptides within all genes
#' @param traces Object of class traces.
#' @param method Which correlation metric to use.
#' @param logtransform Logical, whether to log transform intensities
#' @return List of matrices with pairwise correllations per protein
#' @export
calculateGeneCorrMatrices <- function(traces,
                                  method = c("pearson", "kendall", "spearman"),
                                  logtransform = FALSE){
  genePeptideList <- getGenePepList(traces)

  corrMatrices <- lapply(genePeptideList, function(gene){
    if(logtransform){
      gene <- log1p(gene)
    }
    genecorr <- cor(gene, method = method)
    return(genecorr)
  })
  names(corrMatrices) <- names(genePeptideList)
  traces$geneCorrMatrices <- corrMatrices
  return(traces)
}


#' Filter peptides for having
#' at least one high correlating sibling peptide
#' @param traces Object of class traces.
#' @param cutoff Numeric between 0 and 1. Minimum correlation of a peptide
#' width any sibling peptide.
#' @param filter_by Filtering criterium, maxcorr (default) filters with a cutoff
#' on the correlation, minpval filters with a cutoff on the correlation p-value
#' @param logtransform Logical, whether to log transform intensities
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

filterByMaxCorr <- function(traces, cutoff = 0.85, filter_by = c("maxcorr", "minpval"),
                            logtransform = FALSE, plot = FALSE, PDF=FALSE, name="maxCorrHist", ...) {
  filter_by <- match.arg(filter_by)

  genePeptideList <- getGenePepList(traces)
  maxCorrMatrices <- lapply(genePeptideList, function(gene){

    if(logtransform){
      gene <- log1p(gene)
      gene[gene == -Inf] <- NA # Exclude 0s from the analysis
    }
    genecorr <- cor(gene, use = "pairwise.complete")
    diag(genecorr) <- NA
    maxcorr <- matrixStats::rowMaxs(genecorr, na.rm = T)
    names(maxcorr) <- rownames(genecorr)
    if(filter_by == "maxcorr"){
      return(maxcorr)
    }
    ## Calculate correlation significance
    df = colSums(!is.na(gene)) - 2
    statistic <- sqrt(df) * genecorr/sqrt(1 - genecorr^2)
    pval = 2 * pmin(pt(statistic, df), pt(statistic, df, lower.tail = FALSE))
    pval_adj = matrix(p.adjust(pval), nrow = dim(pval)[1], ncol = dim(pval)[2])
    minp <- matrixStats::rowMins(pval, na.rm = T)
    names(minp) <- rownames(genecorr)
    return(minp)
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
                geom_vline(xintercept = cutoff, color = "red") +
                theme_classic()
    plot(p)
    if (PDF) {
      dev.off()
    }
  }

  if(filter_by == "minpval"){
    filterpeps <- names(allMaxCorrs)[allMaxCorrs < cutoff]
  } else if(filter_by == "maxcorr"){
    filterpeps <- names(allMaxCorrs)[allMaxCorrs > cutoff]
  }

  traces_filt <- subset(traces, filterpeps, "id")
  return(traces_filt)
}

#' @describeIn filterByMaxCorr Filter peptides for having
#' at least one high correlating sibling peptide
#' @export
filterByMaxCorr.tracesList <- function(tracesList, cutoff = 0.85,
                        plot = FALSE, PDF=FALSE, name="maxCorrHist", acrossConditions=TRUE, ...) {

  .tracesListTest(tracesList)
  if (acrossConditions) {
    traces_combined <- combineTracesMutiCond(tracesList)
    traces_maxCorr <- filterByMaxCorr(traces_combined,
                                      cutoff = cutoff,
                                      plot = plot,
                                      PDF=PDF, name=name)
    res <- subset(tracesList, trace_subset_ids = unique(traces_maxCorr$trace_annotation$id))
  } else {
    if (PDF) {pdf(paste0(name,".pdf"),width=3,height=3)}
    res <- lapply(tracesList, filterByMaxCorr.traces,
                  cutoff = cutoff,
                  plot = plot, PDF=F, name=name, ...)
    if (PDF) {dev.off()}
    class(res) <- "tracesList"
  }
  .tracesListTest(res)
  return(res)
}

#' Filter iteratively by maxCorr until all outliers are removed
#' @param traces Object of class traces or tracesList.
#' @param cutoff Numeric between 0 and 1. Minimum correlation of a peptide
#' width any sibling peptide.
#' @param plot logical,wether to print SibPepCorr density plot to R console.
#' Deafult is \code{FALSE}.
#' @param PDF logical, wether to print SibPepCorr density plot to a PDF file.
#' Deafult is \code{FALSE}.
#' @param name Character string with name of the plot, only used if
#' '\code{PDF=TRUE}.PDF file is saved under name.pdf. Default is "maxCorrHist".
#' @return Object of class traces filtered for peptide correlation.
#' @export
iterativeMaxCorrFilter <- function(traces, cutoff = 0.85,
                            plot = FALSE, PDF=FALSE, name="maxCorrHist"){
  UseMethod("iterativeMaxCorrFilter", traces)
}

#' @describeIn iterativeMaxCorrFilter Filter peptides for having
#' at least one high correlating sibling peptide
#' @export
iterativeMaxCorrFilter.traces <- function(traces, cutoff = 0.85,
                        plot = FALSE, PDF=FALSE, name="maxCorrHist", ...) {
  print(nrow(traces$trace_annotation))
  res <- list()
  res[[1]] <- traces
  i=2
  res[[i]] <- filterByMaxCorr(res[[1]], cutoff = cutoff,
                          plot = plot, PDF=PDF, name=name)
  while (!identical(summary(res[[i-1]]),summary(res[[i]]))) {
    print(i)
    print(nrow(res[[i]]$trace_annotation))
    i=i+1
    res[[i]] <- filterByMaxCorr(res[[i-1]], cutoff = cutoff,
                            plot = plot, PDF=PDF, name=name)
  }
  return(res[[i]])
}

#' @describeIn iterativeMaxCorrFilter Filter peptides for having
#' at least one high correlating sibling peptide
#' @export
iterativeMaxCorrFilter.tracesList <- function(traces, cutoff = 0.85,
                        plot = FALSE, PDF=FALSE, name="maxCorrHist", ...) {
  print(nrow(traces[[1]]$trace_annotation))
  print(nrow(traces[[2]]$trace_annotation))
  res <- list()
  res[[1]] <- traces
  i=2
  res[[i]] <- filterByMaxCorr(res[[1]], cutoff = cutoff,
                          plot = plot, PDF=PDF, name=name)
  while (!identical(summary(res[[i-1]]),summary(res[[i]]))) {
    print(i)
    print(nrow(res[[i]][[1]]$trace_annotation))
    print(nrow(res[[i]][[2]]$trace_annotation))
    i=i+1
    res[[i]] <- filterByMaxCorr(res[[i-1]], cutoff = cutoff,
                            plot = plot, PDF=PDF, name=name)
  }
  return(res[[i]])
}

#' Calculate minimum correlation among all peptides that are in a cluster
#' with cluster id > 0
#' @param traces Object of class traces or tracesList.
#' @param plot logical,wether to print SibPepCorr density plot to R console.
#' Deafult is \code{FALSE}.
#' @param PDF logical, wether to print SibPepCorr density plot to a PDF file.
#' Deafult is \code{FALSE}.
#' @param name Character string with name of the plot, only used if
#' '\code{PDF=TRUE}.PDF file is saved under name.pdf. Default is "minCorrHist".
#' @return Object of class traces with calculated minimum.
#' @export
calculateMinCorrClustered.traces <- function(traces, FeatureFinding = NULL,
                                    plot = FALSE, PDF=FALSE, name="minCorrHist", ...) {
  cor_mat <- traces$geneCorrMatrices
  mincorr <- sapply(names(cor_mat), function(prot){
    ann <- traces$trace_annotation[protein_id == prot, .(id, cluster)]
    clustered_peps <- ann[cluster != 0, id]
    mincorr <- min(cor_mat[[prot]][clustered_peps,clustered_peps])
    mincorr
  })
  mincorr[is.infinite(mincorr)] <- NA
  mincorr[mincorr < 0] <- 0
  traces$trace_annotation$mincorr <- mincorr[traces$trace_annotation$protein_id]
  return(traces)
}

#' Calculate minimum correlation of each peptide
#' to all of its sibling peptides
#' @param traces Object of class traces or tracesList.
#' @param plot logical,wether to print SibPepCorr density plot to R console.
#' Deafult is \code{FALSE}.
#' @param PDF logical, wether to print SibPepCorr density plot to a PDF file.
#' Deafult is \code{FALSE}.
#' @param name Character string with name of the plot, only used if
#' '\code{PDF=TRUE}.PDF file is saved under name.pdf. Default is "minCorrHist".
#' @return Object of class traces with calculated minimum.
#' @export
calculateMinCorr <- function(traces, FeatureFinding = NULL,
                            plot = FALSE, PDF=FALSE, name="minCorrHist", ...){
  UseMethod("calculateMinCorr", traces)
}

#' #' @describeIn calculateMinCorr Calculate minimum correlation of each peptide
#' #' to all of its sibling peptides
#' #' @export
calculateMinCorr.traces <- function(traces, FeatureFinding = NULL,
                            plot = FALSE, PDF=FALSE, name="minCorrHist", ...) {
  if (is.null(FeatureFinding)){
    genePeptideList <- getGenePepList(traces)
  } else {
    genePeptideList <- getGenePepListFeatured(traces,FeatureFinding)
  }
  minCorrMatrices <- lapply(genePeptideList, function(gene){
    genecorr <- cor(gene)
    mincorr <- min(genecorr)
    minpeps <- rownames(which(genecorr == min(genecorr), arr.ind=T))
    ints <- apply(gene, 2, max)
    minint <- min(ints)
    if (!is.null(minpeps)){
      data.table(mincorr = mincorr, minint = minint,
                 minpep1 = minpeps[1], minpep2 = minpeps[2])[]
    }
  })
  minCorrMatrices <- do.call(rbind, minCorrMatrices)
  minCorrMatrices[,protein_id := names(genePeptideList)]
  if (plot == TRUE) {
    if (PDF) {
      pdf(paste0(name,".pdf"),width=3,height=3)
    }
    p <- ggplot(data.table(minCorr=minCorrMatrices$mincorr),
                aes(x=minCorr)) +
                geom_histogram(bins=30) +
                theme_classic()
    plot(p)
    #hist(minCorrMatrices$mincorr, breaks = 100)
    if (PDF) {
      dev.off()
    }
  }
  # return(minCorrMatrices[])
  traces$trace_annotation <- merge(traces$trace_annotation,
    minCorrMatrices,by="protein_id",sort=F)
  return(traces)
}

#' #' @describeIn calculateMinCorr Calculate minimum correlation of each peptide
#' #' to all of its sibling peptides
#' #' @export
calculateMinCorr.tracesList <- function(tracesList, FeatureFinding = NULL,
                        plot = FALSE, PDF=FALSE, name="minCorrHist", ...) {

  .tracesListTest(tracesList)
  if (PDF) {pdf(paste0(name,".pdf"))}
  res <- lapply(tracesList, calculateMinCorr.traces,
                plot = plot, PDF=F, name=name, ...)
  if (PDF) {dev.off()}
  class(res) <- "tracesList"
  .tracesListTest(res)
  return(res)
}


#' Fit appropriate distributon to minimum correlation of each peptide
#' to all of its sibling peptides
#' @param traces Object of class traces or tracesList.
#' @param prot_names Vactor with protein names to use for fitting. Default NULL.
#' @param split_value Numeric between 0 and 1. Default is 0.2.
#' @param plot logical,wether to print SibPepCorr density plot to R console.
#' Deafult is \code{FALSE}.
#' @param PDF logical, wether to print SibPepCorr density plot to a PDF file.
#' Deafult is \code{FALSE}.
#' @param name Character string with name of the plot, only used if
#' '\code{PDF=TRUE}.PDF file is saved under name.pdf. Default is "DistrFitted".
#' @return Fitted appropriate distribution. Default is gamma distribution.
#' @export
fitDist <- function(traces, prot_names = NULL,
                      split_value = 0.2, plot = FALSE,
                      PDF=FALSE, name="DistrFitted",...) {
  if (!is.null(prot_names)){
    traces <- subset(traces,
      trace_subset_ids=prot_names,trace_subset_type="protein_id")
  }
  if(! "mincorr" %in% names(traces$trace_annotation)){
    message("no mincorr available: calculating mincorr...")
    traces <- calculateMinCorr(traces, plot = F, PDF=F)
  }
  mincorrs <- as.vector(unique(subset(traces$trace_annotation,
                                      select=c("protein_id","mincorr")))$mincorr)
  DistFitted <- fitdist(1 - mincorrs[mincorrs > split_value], distr=distr)
  if (plot == TRUE) {
    if (PDF) {
      pdf(paste0(name,".pdf"))
    }
    plot(DistFitted)
    if (PDF) {
      dev.off()
    }
  }
  return(DistFitted)
}

#' Estimate p-values of whether a gene is likely to have >1 proteoforms
#' @param traces Object of class traces or tracesList.
#' @param distr Character string which distribution to fit. Default "beta".
#' @param prot_names Vactor with protein names to use for fitting. Default NULL.
#' @param adj.method Character string which p-value adjustment to use.
#' Default "fdr".
#' @param plot logical,wether to print SibPepCorr density plot to R console.
#' Deafult is \code{FALSE}.
#' @param PDF logical, wether to print SibPepCorr density plot to a PDF file.
#' Deafult is \code{FALSE}.
#' @param name Character string with name of the plot, only used if
#' '\code{PDF=TRUE}.PDF file is saved under name.pdf. Default is "DistrFitted".
#' @return Fitted appropriate distribution. Default is gamma distribution.
#' @export
estimateProteoformPval <- function(traces, distr = "beta",
                          prot_names = NULL, adj.method = "fdr",
                          plot = FALSE, PDF=FALSE, name="SplicePval", 
                          design_matrix = NULL, min_sig=2, FDR_cutoff=0.05, ...){
  UseMethod("estimateProteoformPval", traces)
}

#' @describeIn estimateProteoformPval Estimate p-values of whether a
#' gene is likely to have >1 proteoforms
#' @export
estimateProteoformPval.traces <- function(traces, distr = "beta",
                          prot_names = NULL, adj.method = "fdr",
                          plot = FALSE, PDF=FALSE, name="SplicePval", 
                          design_matrix = NULL, min_sig=2, FDR_cutoff=0.05, ...) {
  if (!is.null(prot_names)){
    traces <- subset(traces, trace_subset_ids=prot_names,
      trace_subset_type="protein_id")
  }
  if(! "mincorr" %in% names(traces$trace_annotation)){
    message("no mincorr available: calculating mincorr...")
    traces <- calculateMinCorr(traces, plot = plot, PDF=PDF, name=paste0(name, "_minCorr"))

  } 
  mincorrs <- as.vector(unique(subset(traces$trace_annotation,
                                      select=c("protein_id","mincorr")))$mincorr)
  environment(fitDist) <- environment()
  distr_fitted <- fitDist(traces, plot=plot, PDF=PDF, distr=distr,...)
  # previously we use pgamma directly in the function for our gamma distr 
  # example while we are using a general way to enable utilization of 
  # all ditributions, i.e. "p + distr abbr."
  distrfun <- get(paste("p", distr, sep=""))
  pval <- 1-(distrfun(1 - mincorrs, distr_fitted$estimate[1][[1]],
                   distr_fitted$estimate[2][[1]]))
  pval_adj <- p.adjust(pval, adj.method)
  prot_ids <- unique(traces$trace_annotation$protein_id)
  pval_df <- data.table(prot_ids,pval)
  colnames(pval_df) <- c("protein_id","proteoform_pval")
  pval_adj_df <- data.table(prot_ids,pval_adj)
  colnames(pval_adj_df) <- c("protein_id","proteoform_pval_adj")
  

  if (plot == TRUE) {
    if (PDF) {
      pdf(paste0(name,".pdf"),width=3,height=3)
    }
    p <- ggplot(data.table(p_value=pval),
                aes(x=p_value)) +
                geom_histogram(bins=50) +
                theme_classic()
    plot(p)
    
    q <- ggplot(data.table(adj_p_value=pval_adj),
                aes(x=adj_p_value)) +
                geom_histogram(bins=50) +
                theme_classic()
    plot(q)
    #hist(pval, breaks = 50)
    #hist(pval_adj, breaks = 100)
    if (PDF) {
      dev.off()
    }
  }

  traces$trace_annotation <- merge(traces$trace_annotation, pval_df, by="protein_id", sort=F)
  traces$trace_annotation <- merge(traces$trace_annotation, pval_adj_df, by="protein_id", sort=F)

  return(traces)
}

#' @describeIn estimateProteoformPval Estimate p-values of whether a
#' gene is likely to have >1 proteoforms
#' @export
estimateProteoformPval.tracesList <- function(tracesList, distr = "beta",
                          prot_names = NULL, adj.method = "fdr",
                          plot = FALSE, PDF=FALSE, name="SplicePval", 
                          design_matrix = NULL, min_sig=2, FDR_cutoff=0.05, ...) {

  .tracesListTest(tracesList)
  if (is.null(design_matrix)) {
    stop("Please specifiy a design_matrix.")
  }
  
  proteoform_ann <- data.table("protein_id"=character(),"mincorr"=numeric(),"proteoform_pval"=numeric(),"proteoform_pval_adj"=numeric())
  
  for (r in unique(design_matrix$Replicate)) {
    traces <- tracesList[paste0(unique(design_matrix$Condition),r)]
    class(traces) <- "tracesList"
    traces_combined <- combineTracesMutiCond(traces)
    traces_combined <- filterSinglePeptideHits(traces_combined)
    
    traces_combined_minCorr <- calculateMinCorr(traces_combined,
                                                plot = plot, PDF=PDF,
                                                name=paste0(name, "_minCorr",r))
    
    traces_combined_minCorr_pval <- estimateProteoformPval.traces(traces_combined_minCorr,
                                                                 distr = distr,
                                                                 plot = plot, PDF=PDF,
                                                                 name=paste0(name, "_SplicePval",r))
    
    proteoform_ann <- rbind(proteoform_ann,
                            unique(subset(traces_combined_minCorr_pval$trace_annotation, 
                                          select=c("protein_id","mincorr","proteoform_pval","proteoform_pval_adj"))))
  }
  
  proteoform_ann[, n_prot := .N, by="protein_id"]
  proteoform_ann[, sig := ifelse(proteoform_pval_adj <= FDR_cutoff,1,0)]
  proteoform_ann[, n_sig := sum(sig), by="protein_id"]
  
  proteoform_ann[, medianMinCorr := median(mincorr), by="protein_id"]
  proteoform_ann[, medianPval := median(proteoform_pval), by="protein_id"]
  proteoform_ann[, medianPvaladj := median(proteoform_pval_adj), by="protein_id"]
  
  proteoform_ann[, medianPvaladj := ifelse(n_sig < min_sig, 1, medianPvaladj)]
  
  proteoform_annotation <- unique(subset(proteoform_ann, select=c("protein_id","medianMinCorr","medianPval","medianPvaladj")))
  setnames(proteoform_annotation,c("medianMinCorr","medianPval","medianPvaladj"),c("mincorr","proteoform_pval","proteoform_pval_adj"))
  
  pepTracesList_minCorr_pval <- lapply(tracesList, function(x){
    x$trace_annotation <- merge(x$trace_annotation,
                                proteoform_annotation,by=c("protein_id"),sort=F,all.x=T,all.y=F)
    x
  })
  class(pepTracesList_minCorr_pval) <- "tracesList"

  .tracesListTest(pepTracesList_minCorr_pval)
  return(pepTracesList_minCorr_pval)
}

#' Cluster peptides of gene to unique proteoforms
#' @param traces Object of class traces or tracesList.
#' @param method Character string defining method for clustering
#' Default is "complete".
#' @param plot logical,wether to print SibPepCorr density plot to R console.
#' Deafult is \code{FALSE}.
#' @param PDF logical, wether to print SibPepCorr density plot to a PDF file.
#' Deafult is \code{FALSE}.
#' @param name Character string with name of the plot, only used if
#' '\code{PDF=TRUE}.PDF file is saved under name.pdf. Default is "hclust".
#' @return List including hclust objects for each protein.
#' @export
clusterPeptides <- function(traces,
                            method = c("complete","single","average"), 
                            plot = FALSE, PDF=FALSE, name="hclust", ...) {
  UseMethod("clusterPeptides", traces)
}

#' @describeIn clusterPeptides Cluster peptides of gene to unique proteoforms
#' @export
clusterPeptides.traces <- function(traces,
                             method = c("complete","single","average"),
                             plot = FALSE,
                             PDF = FALSE,
                             name = "hclust",
                             ...){

  method <- match.arg(method)
  # Check if correlations are already computed
  if(! "geneCorrMatrices" %in% names(traces)){
    message("no geneCorrMatrices available: calculating ...")
    traces <- calculateGeneCorrMatrices(traces, ...)
  }
  corMatrices <- traces$geneCorrMatrices

  # Hierarchical clustering
  if (PDF) {
    pdf(paste0(name,".pdf"), width=3, height=3)
  }
  cls <- lapply(1:length(corMatrices), function(i){
    genecorr <- corMatrices[[i]]
    gene_name <- names(corMatrices)[[i]]
    cl <- hclust(as.dist(1-genecorr), method = method, ...)
    if (plot == TRUE) {
      par(cex=0.3, mar=c(18, 5, 5, 1))
      plot(as.dendrogram(cl), ylim = c(0,1), xlab="", ylab="", main="", sub="", axes=FALSE)
      par(cex=1)
      title(xlab="", ylab="", main=gene_name)
      axis(2)
    }
    return(cl)
  })

  if (PDF) {
    dev.off()
  }

  names(cls) <- names(corMatrices)
  traces$peptideClustering <- cls
  return(traces)
}

#' @describeIn clusterPeptides Cluster peptides of gene to unique proteoforms
#' @param traces Object of class traces or tracesList.
#' @return Object of class traces with proteoform annotation
#' @export
clusterPeptides.tracesList <- function(traces, 
                                   method = c("complete","single","average"), 
                                   plot = FALSE, PDF=FALSE, name="hclust", ...) {
  method <- match.arg(method)
  pepTracesList_combi <- combineTracesMutiCond(traces)
  clustering <- clusterPeptides(pepTracesList_combi, method = method,
                                plot = plot, PDF=PDF, 
                                name=name)
  return(clustering)
}

#' Cut clusters from clusterPeptides into N clusters.
#' @param clusterList List of hclust objects generated by clusterPeptides.
#' @param clusterN Number of clusters. Default is 2.
#' @param min_peptides_per_cluster Minimum number of peptides per cluster. Default is 2.
#' @return data.table with information of assigned clusters for each protein and peptide.
#' @export
cutClustersInN <- function(clusterList, clusterN = 2, min_peptides_per_cluster = 2){
  clust <- lapply(seq_along(clusterList), function(idx){
    clusters <- clusterList[[idx]]
    gene_name <- names(clusterList)[[idx]]
    groups <- cutree(clusters, k = clusterN)
    groupsDT <- data.table(id=names(groups),cluster=groups,protein_id=gene_name,
                           proteoform_id = paste0(gene_name,"_",groups))
    groupsDT[, n_pep := .N, by=cluster]
    min_n <- min(groupsDT$n_pep)
    if(min_n < min_peptides_per_cluster){
      groupsDT[, cluster := 1]
      groupsDT[, proteoform_id := paste0(gene_name,"_",1)]
    }
    groupsDT[, n_pep := NULL]
    return(groupsDT)
  })
  
  clustDT <- do.call(rbind, clust)
  
  return(clustDT)
}

#' Cut clusters from clusterPeptides at a specified hight cutoff.
#' @param clusterList List of hclust objects generated by clusterPeptides.
#' @param clusterH Cluster hight cutoff.
#' @return data.table with information of assigned clusters for each protein and peptide.
#' @export
cutClusters <- function(clusterList, clusterH = 0.5){
  clust <- lapply(seq_along(clusterList), function(idx){
    clusters <- clusterList[[idx]]
    gene_name <- names(clusterList)[[idx]]
    groups <- cutree(clusters, h = clusterH)
    groupsDT <- data.table(id=names(groups),cluster=groups,protein_id=gene_name,
                           proteoform_id = paste0(gene_name,"_",groups))
    return(groupsDT)
  })
  
  clustDT <- do.call(rbind, clust)
  
  return(clustDT)
}

#' Annotate traces with cluster information
#' @param traces Object of class traces or tracesList.
#' @param cluster_annotation data.table generated by cutClusters.
#' @param FDR_cutoff FDR cutoff for proteoform significance.Numeric between 0 and 1. 
#' Default is 0.05 (5%).
#' @return Traces with annotated proteoform clusters.
#' @export
annotateTracesWithProteoformClusters <- function(traces, cluster_annotation, FDR_cutoff = 0.05){
  UseMethod("annotateTracesWithProteoformClusters", traces)
}

#' @describeIn annotateTracesWithProteoformClusters Annotate traces with cluster information
#' @export
annotateTracesWithProteoformClusters.traces <- function(traces, cluster_annotation, FDR_cutoff = 0.05){
  traces$trace_annotation <- merge(traces$trace_annotation,cluster_annotation,by=c("id","protein_id"),sort=F,all.x=T,all.y=F)
  traces$trace_annotation[,proteoform_id:=ifelse(is.na(proteoform_id),0,proteoform_id)]
  traces$trace_annotation[,proteoform_id := ifelse(proteoform_pval_adj > FDR_cutoff, protein_id, proteoform_id)]
  traces$trace_annotation[,n_proteoforms:=length(unique(proteoform_id)),by="protein_id"]
  traces$trace_annotation[,n_proteoforms:=ifelse(length(grep("_0",unique(proteoform_id))) == 0, n_proteoforms, n_proteoforms-1),by="protein_id"]
  traces$trace_annotation[,proteoform_id := ifelse(n_proteoforms == 1, protein_id, proteoform_id)]
  .tracesTest(traces)
  return(traces)
}

#' @describeIn annotateTracesWithProteoformClusters Annotate traces with cluster information
#' @export
annotateTracesWithProteoformClusters.tracesList <- function(traces, cluster_annotation, FDR_cutoff = 0.05){
  .tracesListTest(traces)
  combi <- integrateTraceIntensities(traces,design_matrix = NULL,integrate_within = NULL,aggr_fun = "sum")
  combi_ann <- annotateTracesWithProteoformClusters.traces(combi,
                                                           cluster_annotation=cluster_annotation,
                                                           FDR_cutoff=FDR_cutoff)
  
  traces_ann <- lapply(traces, function(x){
    x$trace_annotation <- subset(x$trace_annotation, select = c("protein_id", "id"))
    x$trace_annotation <- merge(x$trace_annotation,
                                combi_ann$trace_annotation,by=c("protein_id", "id"),sort=F,all.x=T,all.y=F)
    x
  })
  class(traces_ann) <- "tracesList"
  
  .tracesListTest(traces_ann)
  return(traces_ann)
}

#' Remove proteins with single peptides
#' @param traces Object of class traces or tracesList.
#' @return Object of class traces containing only proteins
#' with multiple peptides
#' @export
filterSinglePeptideHits <- function(traces){
  UseMethod("filterSinglePeptideHits", traces)
}

#' @describeIn filterSinglePeptideHits Remove proteins with single peptides
#' @export
filterSinglePeptideHits.traces <- function(traces){
  traces$trace_annotation[,n_peptides := .N, by="protein_id"]
  mult_pep_proteins <- unique(traces$trace_annotation[n_peptides >1]$protein_id)
  traces <- subset(traces, trace_subset_ids=mult_pep_proteins,
    trace_subset_type="protein_id")
  return(traces)
}

#' @describeIn filterSinglePeptideHits Remove proteins with single peptides
#' @export
filterSinglePeptideHits.tracesList <- function(tracesList){
  .tracesListTest(tracesList)
  res <- lapply(tracesList, filterSinglePeptideHits.traces)
  class(res) <- "tracesList"
  .tracesListTest(res)
  return(res)
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
    cols <- cols[which(!cols %in% c("SibPepCorr","RepPepCorr","n_peptides", "detectedIn"))]
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
