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

#' Calculates a score for proteoforms based on the difference of within cluster distances
#' and between cluster distances.
#' @param traces Object of class traces with `geneCorrMatrices` and clustered peptides.
#' @param summary_fun the function to use for distance aggregation. Default: mean
#' @return traces object with proteoform_score column in trace_annotation.
#' @export
calculateProteoformScore <- function(traces, summary_fun = "mean"){

  # Check if correlations are already computed
  if(! "geneCorrMatrices" %in% names(traces)){
    stop("no geneCorrMatrices available. Please compute first...")
  }
  # Check if clustering already done
  if(! "cluster" %in% names(traces$trace_annotation)){
    stop("No clustering found please call cutClustersDynamic first.")
  }

  n_fractions <- nrow(traces$fraction_annotation)

  res <- sapply(names(traces$geneCorrMatrices), function(prot){
    clust <- traces$trace_annotation[protein_id == prot]

    if(any(clust$cluster != 100)){
      mat <- traces$geneCorrMatrices[[prot]]
      cl <- unique(clust$cluster)
      cl <- cl[cl != 100]
      if(length(cl) > 2){
        message(paste0("Protein ", prot, " has > 2 clusters."))
        stat_v <- c()
        for (i in seq(1,length(cl),1)){
          cl_i <- cl[cl!=cl[i]]
          cl1 <- clust[clust$cluster == cl_i[1]]$id
          cl2 <- clust[clust$cluster == cl_i[2]]$id
          mat_i <- mat[c(cl1, cl2), c(cl1, cl2)]
          cross <- mat_i[cl1,cl2]
          mat_i[upper.tri(mat_i, diag = T)] <- NA
          ## stat_all <- mean(mat_i, na.rm = T)
          within_c1 <- mat_i[cl1,cl1]
          within_c2 <- mat_i[cl2,cl2]
          stat_within_c1 <- do.call(summary_fun, list(within_c1, na.rm = T))
          stat_within_c2 <- do.call(summary_fun, list(within_c2, na.rm = T))
          stat_within <- min(c(stat_within_c2, stat_within_c1))
          stat_across <- do.call(summary_fun, list(cross, na.rm = T))
          diff_stat <- stat_within - stat_across
          #
          z_stat_within <- atanh(stat_within)
          z_stat_across <- atanh(stat_across)
          z_diff_stat <- z_stat_within - z_stat_across
          dz <- z_diff_stat/(sqrt((1/(n_fractions-3)) + (1/(n_fractions-3))))
          pval <- 2*(1 - pnorm(abs(dz)))
          #
          stat_v <- append(stat_v,list(list(diff_stat, z_diff_stat, dz, pval)))
        }
        sel_min_diff <- which(unlist(lapply(stat_v, function(l){l[1]})) == min(unlist(lapply(stat_v, function(l){l[1]})), na.rm = T))[1]
        return(unlist(stat_v[[sel_min_diff]]))
      } else {
        cl1 <- clust[clust$cluster == cl[1]]$id
        cl2 <- clust[clust$cluster == cl[2]]$id
        mat <- mat[c(cl1, cl2), c(cl1, cl2)]
        cross <- mat[cl1,cl2]
        mat[upper.tri(mat, diag = T)] <- NA
        ## stat_all <- mean(mat, na.rm = T)
        within_c1 <- mat[cl1,cl1]
        within_c2 <- mat[cl2,cl2]
        stat_within_c1 <- do.call(summary_fun, list(within_c1, na.rm = T))
        stat_within_c2 <- do.call(summary_fun, list(within_c2, na.rm = T))
        stat_within <- min(c(stat_within_c2, stat_within_c1))
        stat_across <- do.call(summary_fun, list(cross, na.rm = T))
        diff_stat <- stat_within - stat_across
        #
        z_stat_within <- atanh(stat_within)
        z_stat_across <- atanh(stat_across)
        z_diff_stat <- z_stat_within - z_stat_across
        dz <- z_diff_stat/(sqrt((1/(n_fractions-3)) + (1/(n_fractions-3))))
        pval <- 2*(1 - pnorm(abs(dz)))
        #
        return(list(diff_stat, z_diff_stat, dz, pval))
      }
    } else{
      return(NA)
    }
  })
  ## Add to trace annotation
  #traces$trace_annotation$proteoform_score <- res[traces$trace_annotation$protein_id]
  traces$trace_annotation$proteoform_score <- unlist(lapply(res[traces$trace_annotation$protein_id], function(l){l[1]}))
  traces$trace_annotation$proteoform_score_z <- unlist(lapply(res[traces$trace_annotation$protein_id], function(l){l[2]}))
  traces$trace_annotation$proteoform_score_dz <- unlist(lapply(res[traces$trace_annotation$protein_id], function(l){l[3]}))
  traces$trace_annotation$proteoform_score_pval <- unlist(lapply(res[traces$trace_annotation$protein_id], function(l){l[4]}))
  pval_dt <- unique(traces$trace_annotation[,.(proteoform_score_pval), by=.(protein_id)])
  pval_dt[, proteoform_score_pval_adj := p.adjust(proteoform_score_pval, method = "BH")]
  traces$trace_annotation <- merge(traces$trace_annotation, pval_dt, by=c("protein_id", "proteoform_score_pval"), sort = F)
  return(traces)
}

#' Plots p-value histograms for the proteoform p-values.
#' @param traces Object of class traces with proteoform scores and p-values.
#' @param name Name of the generated PDF.
#' @param PDF logical, wether to print plot to a PDF file.
#' @return plot
#' @export
plotProteoformPvalHist <- function(traces, name="proteoform_pval_histograms", PDF=TRUE){
  if (PDF){
    pdf(paste0(name,".pdf"), width=4, height=4)
  }
    scores <- unique(traces$trace_annotation[,.(proteoform_score,proteoform_score_pval,proteoform_score_pval_adj), by=.(protein_id)])
    scores <- scores[!is.na(scores$proteoform_score)]
    p <- ggplot(scores, aes(x=proteoform_score_pval)) +
      geom_histogram(bins = 100) +
      theme_classic()
    print(p)
    p <- ggplot(scores, aes(x=proteoform_score_pval_adj)) +
      geom_histogram(bins = 80) +
      theme_classic()
    print(p)
  if (PDF){
    dev.off()
  }
}

#' Plots pseudo volcano for the proteoform scores and p-values.
#' @param traces Object of class traces with proteoform scores and p-values.
#' @param name Name of the generated PDF.
#' @param PDF logical, wether to print plot to a PDF file.
#' @return plot
#' @export
plotProteoformVolcano <- function(traces, name="proteoform_pseudo_volcano", PDF=TRUE){
  if (PDF){
    pdf(paste0(name,".pdf"), width=4, height=4)
  }
  scores <- unique(traces$trace_annotation[,.(proteoform_score,proteoform_score_pval,proteoform_score_pval_adj), by=.(protein_id)])
  scores <- scores[!is.na(scores$proteoform_score)]
  p <- ggplot(scores, aes(x=proteoform_score, y=-log10(proteoform_score_pval_adj))) +
    geom_point(colour="#003366", alpha=0.5) +
    theme_classic()
  print(p)
  if (PDF){
    dev.off()
  }
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

#' Remove any single peptides that are on the outside of the tree
#' @param tree hclust object
#' @param clusters vector of cluster assignments (as produced by cutree)
#' @param pruning_cutoff int, height up to which single peptides are called outliers
prune_tree <- function(tree, clusters, pruning_cutoff = 0.3){
  last_row <- dim(tree$merge)[1]
  hheight <- tree$height[last_row]
  last_merge <- tree$merge[last_row,]
  if(all(last_merge > 0) | hheight < pruning_cutoff){
    return(clusters)
  }else{
    while(hheight > pruning_cutoff & any(last_merge < 0) & last_row > 0){
      merged_w_single <- last_merge < 0
      del_idx <- -last_merge[merged_w_single]
      clusters[del_idx] <- 0
      last_row <- last_row -1
      if(last_row == 0){break}
      hheight <- tree$height[last_row]
      last_merge <- tree$merge[last_row,]
    }
    return(clusters)
  }
}

#' Compute the silhouette score of every peptide
#' @param tree hclust object
#' @param dist distances that the tree is based on
#' @param clusters vector of cluster assignments (as produced by cutree)
#' @param h int, height at which to cut the tree if clusters are not provided
computeSilhouettes <- function(tree, dist, clusters = NULL, h = 0.3){
  if(is.null(clusters)){
    clusters <- cutree(tree, h = h)
  }
  sh <- cluster::silhouette(x = clusters, dmatrix = 1-dist)
  if(is.na(sh)){
    cluster_info <- data.table(id = names(clusters),
                               cluster = clusters,
                               neighbor =NA,
                               sil_width=NA,
                               cluster_sh = NA)
  }else{
    cluster_info <- data.table(id = names(clusters), sh[,1:3])
    cluster_info[, cluster_sh := mean(sil_width), by=.(cluster)]
  }
  return(cluster_info)
}


#' Compute the height of the tree without the outlier peptides
#' @param tree hclust object
#' @param clusters vector of cluster assignments (as produced by cutree)
get_tree_height <- function(tree, clusters){
  outliers <- which(clusters == 0)
  outlier_ind <- c()
  for (outlier in outliers){
    outlier_ind <- c(which(tree$merge == -outlier, arr.ind = T)[1], outlier_ind)
    while(any(tree$merge == outlier_ind[1])){
      outlier_ind = c(which(tree$merge == outlier_ind[1], arr.ind = T)[1], outlier_ind)
    }
  }
  if(!is.null(outlier_ind)){
    th <- tree$height[-outlier_ind]
  }else{
    th <- tree$height
  }
  return(max(th))
}

optimize_silhouette <- function(tree, dist, k = 2:5){
  means <- numeric(0)
  shs <- list()
  res <- for(x in k){
    clusters <- cutree(tree, k = x)
    sh <- cluster::silhouette(x = clusters, dist = as.dist(dist))
    sw <- summary(sh)$clus.avg.widths
    ## sw <- sw[sw >0] # Don't take single peptides into account
    shs[[x]] <- sh
    means <- c(means, mean(sw))
  }
  n_optim <- which.max(means) + k[1] -1
  return(list(n_optim = n_optim, avg_sil = max(means), silhoutette = shs[[n_optim]]))
}

#' Remove outlier peptides and take the silhouette width into account
#' @param traces Object of class traces or tracesList.
#' @param initial_cut_height int, height at which tree is initially cut
#' @param pruning_cutoff int, height up to which single peptides are called outliers
#' @param PDF logical, wether to print SibPepCorr density plot to a PDF file.
#' Deafult is \code{FALSE}.
#' @param name Character string with name of the plot, only used if
#' '\code{PDF=TRUE}.PDF file is saved under name.pdf. Default is "robust_clustering".
#' @return tracesList including cluster information in trace_annotation
#' @export
robust_clustering <- function(traces, pruning_cutoff=0.3,
                              initial_cut_height= 0.3,
                              method = c("complete","single","average"),
                              PDF=F, name = 'robust_clustering'){
  method <- match.arg(method)
  if (PDF) {
    pdf(paste0(name,".pdf"), width=3, height=3)
  }
  res <- lapply(names(traces$peptideClustering), function(prot){
    tree <- traces$peptideClustering[[prot]]
    dist <- 1-traces$geneCorrMatrices[[prot]]
    initial_clustering <- cutree(tree, h = initial_cut_height)
    #print(prot)
    t <- table(initial_clustering)
    singlepep_clusters <- names(t[t==1])
    ## pruned_clusters <- prune_tree(tree,
    ##                               initial_clustering,
    ##                               pruning_cutoff = pruning_cutoff)
    pruned_clusters <- initial_clustering
    pruned_clusters[pruned_clusters %in% as.numeric(singlepep_clusters)] <- 0
    ## tree_height_adj <- get_tree_height(tree, pruned_clusters)

    remaining_peps <- names(pruned_clusters[pruned_clusters != 0])
    n_peptides <- length(remaining_peps)
    if(n_peptides > 2){
      n_initial <- length(unique(pruned_clusters[pruned_clusters != 0]))
      d_sub <- dist[remaining_peps,remaining_peps]
      tree_sub <- hclust(as.dist(d_sub),
                         method = method)

      optim_max <- min(n_initial+1, n_peptides - 1)
      cluster_sil <- optimize_silhouette(tree_sub, d_sub, k = 2:optim_max)
      cluster_info <- data.table(id = tree_sub$labels,
                                 cluster_sil$silhoutette[,1:3])
      cluster_info[, cluster_sh := mean(sil_width), by=.(cluster)]
      tree_height_adj <- get_tree_height(tree_sub, cluster_info$cluster)
      ## Add back pruned peptides
      outlier_peps <- names(pruned_clusters[pruned_clusters == 0])
      cluster_info_removed <- data.table(id = outlier_peps,
                                         cluster = 0,
                                         neighbor =NA,
                                         sil_width=NA,
                                         cluster_sh = NA)
      cluster_info <- rbind(cluster_info, cluster_info_removed)
      cluster_info[, mean_sh := cluster_sil$avg_sil]
    }else{
      cluster_info <- data.table(id = names(pruned_clusters),
                                 cluster = pruned_clusters,
                                 neighbor =NA,
                                 sil_width=NA,
                                 cluster_sh = NA,
                                 mean_sh = NA)
      tree_height_adj <- NA
    }
    cluster_info$tree_height_adj <- tree_height_adj
    cluster_info <- cluster_info[order(match(id, names(initial_clustering)))]

    ## cluster_info <- computeSilhouettes(tree, dist, clusters = pruned_clusters)

    # Remove additional single peptide clusters
    t <- table(cluster_info$cluster)
    singlepep_clusters <- names(t[t==1])
    cluster_info[cluster %in% as.numeric(singlepep_clusters), cluster := 0]
    cluster_info[,n_clusters := length(unique(cluster[cluster > 0]))]
    ## print(prot)

    if (PDF) {
      par(cex=0.3, mar=c(18, 5, 5, 1))
      tree$labels <- paste(cluster_info$id, cluster_info$cluster, sep="_") # order has to be conserved!
      plot(as.dendrogram(tree), ylim = c(0,1),
           xlab="", ylab="", main="", sub="",
           axes=FALSE)
      ## rect.hclust(tree, cluster = initial_clustering, h = 0.3)
      par(cex=1)
      title(xlab="", ylab="", main=prot)
      axis(2)
    }
    return(cluster_info)
  })
  if (PDF) {
    dev.off()
  }
  cluster_info <- rbindlist(res)
  cluster_info$cluster <- as.factor(cluster_info$cluster)
  traces$trace_annotation <- merge(traces$trace_annotation, cluster_info, by="id", all.x = T, sort = F)
  return(traces)
}

#' Fit appropriate distributon to minimum correlation of each peptide
#' to all of its sibling peptides
#' @param x Vector with values to fit distribution to.
#' @param split_value Numeric between 0 and 1. Default is 0.5.
#' @param distribution String, which distribution to use.
#' @param plot logical,wether to print SibPepCorr density plot to R console.
#' Deafult is \code{FALSE}.
#' @param PDF logical, wether to print SibPepCorr density plot to a PDF file.
#' Deafult is \code{FALSE}.
#' @param name Character string with name of the plot, only used if
#' '\code{PDF=TRUE}.PDF file is saved under name.pdf. Default is "DistrFitted".
#' @return Fitted appropriate distribution. Default is gamma distribution.
#' @export
fitDistCensored <- function(x,
                      split_value = 0.5,
                      start = NULL,
                      distribution = c("beta", "gamma", "f"),
                      plot = FALSE,
                      PDF=FALSE, name="DistrFitted",...) {

  distribution <- match.arg(distribution)

  mincd <- data.frame(left=x, right=x)
  ## print(mincd)
  mincd[mincd$left > split_value,"left"] <- NA

  fit <- fitdistrplus::fitdistcens(censdata = as.data.frame(mincd),
                                   distr = distribution, start = start)
  if (plot == TRUE) {
    if (PDF) {
      pdf(paste0(name,".pdf"))
    }
    estimates <- split(unname(fit$estimate), names(fit$estimate))
    distrdensfun <- paste0('d', fit$distname)
    p <- ggplot(mincd, aes(x=right)) +
      geom_histogram(alpha = 0.5, aes(y=..density..), bins = 50) +
      stat_function(fun = distrdensfun, args = estimates) +
      theme_minimal()
    print(p)
    if (PDF) {
      dev.off()
    }
  }
  return(fit)
}

#' Fit a censored distribution to the minimum correlation of each gene to estimate
#' a pvalue of isoforms.
#' @param traces Object of class traces or tracesList.
#' @param distribution String, which distribution to use.
#' @param prot_names Vactor with protein names to use for fitting. Default NULL.
#' @param adj.method Character string which p-value adjustment to use.
#' Default "fdr".
#' @param split_value Numeric between 0 and 1, where to censor distribution.
#' Default is 0.5.
#' @param plot logical,wether to print SibPepCorr density plot to R console.
#' Deafult is \code{FALSE}.
#' @param PDF logical, wether to print SibPepCorr density plot to a PDF file.
#' Deafult is \code{FALSE}.
#' @param name Character string with name of the plot, only used if
#' '\code{PDF=TRUE}.PDF file is saved under name.pdf. Default is "DistrFitted".
#' @return Fitted appropriate distribution. Default is gamma distribution.
#' @export
estimateProteoformPvalCens <- function(traces, distr = "beta",
                                   prot_names = NULL, adj.method = "fdr",
                                   plot = FALSE, PDF=FALSE, name="SplicePval",
                                   design_matrix = NULL, min_sig=2, FDR_cutoff=0.05, ...){
  UseMethod("estimateProteoformPvalCens", traces)
}

#' @describeIn estimateProteoformPvalCens Estimate p-values of whether a
#' gene is likely to have >1 proteoforms
#' @export
estimateProteoformPvalCens.traces <- function(traces, distribution = c("beta", "gamma"),
                                          prot_names = NULL, adj.method = "fdr",
                                          split_value = 0.5,
                                          plot = FALSE, PDF=FALSE, name="SplicePval",
                                          design_matrix = NULL, min_sig=2, FDR_cutoff=0.05, ...) {

  distribution <- match.arg(distribution)
  if (!is.null(prot_names)){
    traces <- subset(traces,
                     trace_subset_ids=prot_names,trace_subset_type="protein_id")
  }
  if(! "mincorr" %in% names(traces$trace_annotation)){
    stop("no mincorr available")
  }
  minc_table <- traces$trace_annotation[,unique(mincorr), by = "protein_id"]
  minc <- minc_table$V1
  names(minc) <- minc_table$protein_id

  x <- 1-minc
  fit <- fitDistCensored(x[!is.na(x)], distribution = distribution,
                         plot = plot, PDF = PDF,
                         split_value = split_value)

  distrpfun <- paste0('p', fit$distname)
  args <- split(unname(fit$estimate), names(fit$estimate))
  args$q <- x
  args$lower.tail <- FALSE
  pval <- do.call(distrpfun, args)
  names(pval) <- names(minc)
  pval_adj <- p.adjust(pval, method = "fdr")
  names(pval_adj) <- names(minc)

  traces$trace_annotation$pval <- pval[traces$trace_annotation$protein_id]
  traces$trace_annotation$pval_adj <- pval_adj[traces$trace_annotation$protein_id]

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
    if (PDF) {
      dev.off()
    }
  }

    return(traces)
}

#' Cut clusters from clusterPeptides into N clusters with more than 1 peptide.
#' @param traces Object of class traces with a member `peptideClustering` as produced by `clusterPeptides`.
#' @param clusterN How many clusters to cut the tree into. Default: 2.
#' @param min_peptides_per_cluster How many peptides per cluster.
#' @details Cuts the hierarchical tree so that there are N clusters with at least the specified number of peptides.
#' If there are single peptides above the cutoff height they will get the label 100.
#' @return traces object with cluster annotations in trace_annotation
#' @export

cutClustersInNreal <- function(traces, clusterN = 2, min_peptides_per_cluster = 2){

  if("peptideClustering" %in% names(traces)){
    clusterList <- traces$peptideClustering
  }else{
    stop("No peptide clustering found. Please run clusterPeptides first.")
  }

  clust <- lapply(seq_along(clusterList), function(idx){
    tree <- clusterList[[idx]]
    n_peptides <- length(tree$labels)
    n_real_clusters <- 0
    k <- clusterN
    ## print(idx)
    while(n_real_clusters < clusterN){
      groups <- cutree(tree, k = k)
      n_per_cluster <- table(groups)
      is_multipep <- n_per_cluster >= min_peptides_per_cluster
      singlepep_clusters <- as.numeric(names(n_per_cluster[!is_multipep]))
      groups[groups %in% singlepep_clusters] <- 100
      n_real_clusters <- sum(is_multipep)
      k <- k + 1
      if(k >= n_peptides){
        groups[] <- 100
        break
      }
    }
    ## Renumber groups so that they are 1,2,.. 100
    cn <- max(groups)
    clf <- as.factor(groups)
    cll <- levels(clf)
    hasnoisecl <- any(cll == 100)
    cln <- length(cll)
    if(cn != cln){
      for (i in 1:cln) groups[clf == cll[i]] <- i
    }
    if(hasnoisecl) groups[groups == max(groups)] <- 100

    return(groups)
  })

  names(clust) <- names(clusterList)
  #return(clust)

  traces$trace_annotation[, cluster := clust[[protein_id]][id], by=c('protein_id','id')]
  return(traces)
}

#' Cut clusters from clusterPeptides dynamically into clusters with more than 1 peptide.
#' @param traces Object of class traces with a member `peptideClustering` as produced by `clusterPeptides`.
#' @param min_peptides_per_cluster How many peptides per cluster.
#' @details Cuts the hierarchical tree so that there are N clusters with at least the specified number of peptides.
#' The cutreeDynamic function from dynamicTreeCut is applied.
#' If there are single peptides above the cutoff height they will get the label 100.
#' @return traces object with cluster annotations in trace_annotation
#' @export

cutClustersDynamic <- function(traces, min_peptides_per_cluster = 2){

  if("peptideClustering" %in% names(traces)){
    clusterList <- traces$peptideClustering
  }else{
    stop("No peptide clustering found. Please run clusterPeptides first.")
  }

  clust <- lapply(seq_along(clusterList), function(idx){
    tree <- clusterList[[idx]]
    corr <- 1-traces$geneCorrMatrices[[idx]]
    peptides <- tree$labels
    groups <- cutreeDynamic(tree, minClusterSize = 2, distM = corr, cutHeight=0.99, minSplitHeight=0.9, deepSplit=0, pamStage = F)
    names(groups) <- peptides
    return(groups)
  })

  names(clust) <- names(clusterList)

  traces$trace_annotation[, cluster := clust[[protein_id]][id], by=c('protein_id','id')]
  traces$trace_annotation[, cluster := ifelse(cluster==0, 100, cluster)]
  return(traces)
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
#' @param FDR_cutoff FDR cutoff for proteoform significance. Numeric between 0 and 1.
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


#' Annotate traces with a proteoform_id based on
#' a selected score cutoff.
#' @param traces Object of class traces.
#' @param score_cutoff Proteoform score cutoff. Numeric between 0 and 1. Default is 0.1.
#' @param adj_pval_cutoff adjusted p-value cutoff. Numeric between 0 and 1. Default is 0.1.
#' @return Object of class traces
#' @export
annotateTracesWithProteoforms <- function(traces, score_cutoff = 0.1, adj_pval_cutoff = 0.1){
  traces$trace_annotation[,proteoform_id := ifelse(((proteoform_score >= score_cutoff) & (proteoform_score_pval_adj <= adj_pval_cutoff)), paste(protein_id, cluster, sep='_'), protein_id)]
  traces$trace_annotation[,proteoform_id := ifelse(cluster == 100, paste(protein_id, 0, sep='_'), proteoform_id)]
  traces$trace_annotation[,n_proteoforms:=length(unique(proteoform_id)),by="protein_id"]
  traces$trace_annotation[,n_proteoforms:=ifelse(length(grep("_0",unique(proteoform_id))) == 0, n_proteoforms, n_proteoforms-1),by="protein_id"]
  traces$trace_annotation[,proteoform_id := ifelse(n_proteoforms %in% c(0,1), protein_id, proteoform_id)]
  .tracesTest(traces)
  return(traces)
}
