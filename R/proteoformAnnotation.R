#' Filter peptides for having
#' at least one high correlating sibling peptide
#' @param traces Object of class traces or tracesList.
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
filterByMaxCorr <- function(traces, cutoff = 0.85,
                            plot = FALSE, PDF=FALSE, name="maxCorrHist"){
  UseMethod("filterByMaxCorr", traces)
}

#' @describeIn filterByMaxCorr Filter peptides for having
#' at least one high correlating sibling peptide
#' @export
filterByMaxCorr.traces <- function(traces, cutoff = 0.85,
                            plot = FALSE, PDF=FALSE, name="maxCorrHist", ...) {
  genePeptideList <- getGenePepList(traces)
  maxCorrMatrices <- lapply(genePeptideList, function(gene){
    genecorr <- cor(gene)
    genecorr[genecorr == 1] <- NA
    maxcorr <- rowMaxs(genecorr, na.rm = T)
    names(maxcorr) <- rownames(genecorr)
    return(maxcorr)
  })

  names(maxCorrMatrices) <- names(genePeptideList)
  allMaxCorrs <- unlist(maxCorrMatrices)
  names(allMaxCorrs) <- gsub(".*\\.", "", names(allMaxCorrs))
  if (plot == TRUE){
    if (PDF) {
      pdf(paste0(name,".pdf"))
    }
    hist(allMaxCorrs, breaks = 30)
    plot(density(allMaxCorrs))
    if (PDF) {
      dev.off()
    }
  }

  filterpeps <- names(allMaxCorrs)[allMaxCorrs > cutoff]

  traces_filt <- subset(traces, filterpeps, "id")
  return(traces_filt)
}

#' @describeIn filterByMaxCorr Filter peptides for having
#' at least one high correlating sibling peptide
#' @export
filterByMaxCorr.tracesList <- function(tracesList, cutoff = 0.85,
                        plot = FALSE, PDF=FALSE, name="maxCorrHist", ...) {

  .tracesListTest(tracesList)
  if (PDF) {pdf(paste0(name,".pdf"))}
  res <- lapply(tracesList, filterByMaxCorr.traces,
                cutoff = cutoff,
                plot = plot, PDF=F, name=name, ...)
  if (PDF) {dev.off()}
  class(res) <- "tracesList"
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
calculateMinCorr <- function(traces,
                            plot = FALSE, PDF=FALSE, name="minCorrHist", ...){
  UseMethod("calculateMinCorr", traces)
}

#' @describeIn calculateMinCorr Calculate minimum correlation of each peptide
#' to all of its sibling peptides
#' @export
calculateMinCorr.traces <- function(traces,
                            plot = FALSE, PDF=FALSE, name="minCorrHist", ...) {
  genePeptideList <- getGenePepList(traces)
  minCorrMatrices <- lapply(genePeptideList, function(gene){
    genecorr <- cor(gene)
    mincorr <- min(genecorr)
    minpeps <- rownames(which(genecorr == min(genecorr), arr.ind=T))
    ints <- apply(gene, 2, max)
    minint <- min(ints)
    data.table(mincorr = mincorr, minint = minint,
                      minpep1 = minpeps[1], minpep2 = minpeps[2])[]
  })
  minCorrMatrices <- do.call(rbind, minCorrMatrices)
  minCorrMatrices[,protein_id := names(genePeptideList)]
  if (plot == TRUE) {
    if (PDF) {
      pdf(paste0(name,".pdf"))
    }
    hist(minCorrMatrices$mincorr, breaks = 100)
    if (PDF) {
      dev.off()
    }
  }
  # return(minCorrMatrices[])
  traces$trace_annotation <- merge(traces$trace_annotation,
    minCorrMatrices,by="protein_id",sort=F)
  return(traces)
}

#' @describeIn calculateMinCorr Calculate minimum correlation of each peptide
#' to all of its sibling peptides
#' @export
calculateMinCorr.tracesList <- function(tracesList,
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


#' Fit gamma distributon to minimum correlation of each peptide
#' to all of its sibling peptides
#' @param traces Object of class traces or tracesList.
#' @param prot_names Vactor with protein names to use for fitting. Default NULL.
#' @param distr Character string which distribution to fit. Default "gamma".
#' @param split_value Numeric between 0 and 1. Default is 0.2.
#' @param plot logical,wether to print SibPepCorr density plot to R console.
#' Deafult is \code{FALSE}.
#' @param PDF logical, wether to print SibPepCorr density plot to a PDF file.
#' Deafult is \code{FALSE}.
#' @param name Character string with name of the plot, only used if
#' '\code{PDF=TRUE}.PDF file is saved under name.pdf. Default is "GammaFitted".
#' @return Fitted gamma distribution
#' @export
fitGammaDist <- function(traces, prot_names = NULL, distr = "gamma",
                      split_value = 0.2, plot = FALSE,
                      PDF=FALSE, name="GammaFitted") {
  if (!is.null(prot_names)){
    traces <- subset(traces,
      trace_subset_ids=prot_names,trace_subset_type="protein_id")
  }
  if(! "mincorr" %in% names(traces$trace_annotation)){
    message("no mincorr available: calculating mincorr...")
    traces <- calculateMinCorr(traces, plot = F, PDF=F)
  }
  mincorrs <- traces$trace_annotation$mincorr
  GammaFitted <- fitdist(1 - mincorrs[mincorrs > split_value], distr)
  if (plot == TRUE) {
    if (PDF) {
      pdf(paste0(name,".pdf"))
    }
    plot(GammaFitted)
    if (PDF) {
      dev.off()
    }
  }
  return(GammaFitted)
}

#' Estimate p-values of whether a gene is likely to have >1 proteoforms
#' @param traces Object of class traces or tracesList.
#' @param prot_names Vactor with protein names to use for fitting. Default NULL.
#' @param adj.method Character string which p-value adjustment to use.
#' Default "fdr".
#' @param plot logical,wether to print SibPepCorr density plot to R console.
#' Deafult is \code{FALSE}.
#' @param PDF logical, wether to print SibPepCorr density plot to a PDF file.
#' Deafult is \code{FALSE}.
#' @param name Character string with name of the plot, only used if
#' '\code{PDF=TRUE}.PDF file is saved under name.pdf. Default is "GammaFitted".
#' @return Fitted gamma distribution
#' @export
estimateProteoformPval <- function(traces,
                          prot_names = NULL, adj.method = "fdr",
                          plot = FALSE, PDF=FALSE, name="SplicePval", ...){
  UseMethod("estimateProteoformPval", traces)
}

#' @describeIn estimateProteoformPval Estimate p-values of whether a
#' gene is likely to have >1 proteoforms
#' @export
estimateProteoformPval.traces <- function(traces,
                          prot_names = NULL, adj.method = "fdr",
                          plot = FALSE, PDF=FALSE, name="SplicePval", ...) {
  if (!is.null(prot_names)){
    traces <- subset(traces, trace_subset_ids=prot_names,
      trace_subset_type="protein_id")
  }
  if(! "mincorr" %in% names(traces$trace_annotation)){
    message("no mincorr available: calculating mincorr...")
    traces <- calculateMinCorr(traces, plot = F, PDF=F)
  }
  mincorrs <- traces$trace_annotation$mincorr
  distr_fitted <- fitGammaDist(traces, plot=plot, PDF=PDF)
  # here we use pgamma directly in the function for our gamma distr example,
  # if we want to use other distr, we need a general way to represent
  # all ditribution, e.g., "p + distr abbr."
  pval <- 1-(pgamma(1 - mincorrs, shape = distr_fitted$estimate[1][[1]],
                    rate = distr_fitted$estimate[2][[1]]))
  pval_adj <- p.adjust(pval, adj.method)
  if (plot == TRUE) {
    if (PDF) {
      pdf(paste0(name,".pdf"))
    }
    hist(pval, breaks = 50)
    hist(pval_adj, breaks = 100)
    if (PDF) {
      dev.off()
    }
  }
  traces$trace_annotation[, proteoform_pval := pval_adj]
  traces$trace_annotation[, proteoform_pval_adj := pval_adj]
  return(traces)
}

#' @describeIn estimateProteoformPval Estimate p-values of whether a
#' gene is likely to have >1 proteoforms
#' @export
estimateProteoformPval.tracesList <- function(tracesList,
                          prot_names = NULL, adj.method = "fdr",
                          plot = FALSE, PDF=FALSE, name="SplicePval", ...) {

  .tracesListTest(tracesList)
  if (PDF) {pdf(paste0(name,".pdf"))}
  res <- lapply(tracesList, estimateProteoformPval.traces,
                prot_names = prot_names, adj.method = adj.method,
                plot = plot, PDF=F, name=name, ...)
  if (PDF) {dev.off()}
  class(res) <- "tracesList"
  .tracesListTest(res)
  return(res)
}

#' Cluster peptides of gene to unique proteoforms
#' @param traces Object of class traces or tracesList.
#' @param method Character string defining method for clustering
#' Default is "complete".
#' @param clusterH Numeric Cluster hight for cutree. Default is 0.5.
#' @param clusterN Integer Number of clusters. Default is NULL.
#' @param index Character string for method of cluster number estimation
#' Default is "silhouette".
#' @param plot logical,wether to print SibPepCorr density plot to R console.
#' Deafult is \code{FALSE}.
#' @param PDF logical, wether to print SibPepCorr density plot to a PDF file.
#' Deafult is \code{FALSE}.
#' @param name Character string with name of the plot, only used if
#' '\code{PDF=TRUE}.PDF file is saved under name.pdf. Default is "hclust".
#' @return Object of class traces with proteoform annotation
#' @export
clusterPeptides <- function(traces,
                    method = c("complete","single"), clusterH = 0.5,
                    clusterN = NULL, index = "silhouette",
                    plot = FALSE, PDF=FALSE, name="hclust", ...) {
  UseMethod("clusterPeptides", traces)
}

#' @describeIn clusterPeptides Cluster peptides of gene to unique proteoforms
#' @export
clusterPeptides.traces <- function(traces,
                    method = c("complete","single"), clusterH = 0.5,
                    clusterN = NULL, index = "silhouette",
                    plot = FALSE, PDF=FALSE, name="hclust", ...) {
  method <- match.arg(method)
  if (PDF) {
    pdf(paste0(name,".pdf"))
  }
  genePeptideList <- getGenePepList(traces)
  idx <- seq_along(genePeptideList)
  clustMatricesMethod <- lapply(idx, function(i){
    gene <- genePeptideList[[i]]
    gene_name <- names(genePeptideList)[[i]]
    genecorr <- cor(gene)
    if (is.null(clusterH)) {
      if (is.null(clusterN)) {
        nbClust_res <- NbClust(data=genecorr, diss=as.dist(1-genecorr),
        distance = NULL, min.nc=2,
        max.nc=min(4,nrow(genecorr)-1)[1],method = method, index = index)
        clusterN <- nbClust_res$Best.nc[1]
      }
    }
    cl <- hclust(as.dist(1-genecorr), method)
    if (plot == TRUE) {
      #plot(cl, labels = FALSE)
      plot(as.dendrogram(cl), ylim = c(0,1))
    }
    groups <- cutree(cl, k = clusterN, h = clusterH)
    groupsDT <- data.table(id=names(groups),cluster=groups,protein_id=gene_name,
                          proteoform_id = paste0(gene_name,"_",groups))
    return(groupsDT)
  })
  if (PDF) {
    dev.off()
  }
  #return(clustMatricesMethod)
  clustMatrices <- do.call(rbind, clustMatricesMethod)
  traces$trace_annotation <- merge(traces$trace_annotation,
    clustMatrices,by=c("id","protein_id"),sort=F)
  return(traces)
}

#' @describeIn clusterPeptides Cluster peptides of gene to unique proteoforms
#' @param traces Object of class traces or tracesList.
#' @return Object of class traces with proteoform annotation
#' @export
clusterPeptides.tracesList <- function(tracesList,
                    method = c("complete","single"), clusterH = 0.5,
                    clusterN = NULL, index = "silhouette",
                    plot = FALSE, PDF=FALSE, name="hclust", ...) {

  .tracesListTest(tracesList)
  if (PDF) {pdf(paste0(name,".pdf"))}
  res <- lapply(tracesList, clusterPeptides.traces,
                        method = method, clusterH = clusterH,
                        clusterN = clusterN, index = index,
                        plot = plot, PDF=F, name=name, ...)
  if (PDF) {dev.off()}
  class(res) <- "tracesList"
  .tracesListTest(res)
  return(res)
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
  traces_all <- dcast(traces_combi,id ~ variable+cond,value.var="value",fill=0)
  names(traces_all) <- c("id",seq(1,ncol(traces_all)-1,1))
  setcolorder(traces_all, c(seq(1,ncol(traces_all)-1,1),"id"))

  trace_type_all <- tracesList[[1]]$trace_type

  trace_annotation <- lapply(tracesList, function(t){
    cols <- names(tracesList$minus$trace_annotation)
    cols <- cols[which(!cols %in% c("SibPepCorr","RepPepCorr"))]
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
  class(combi_traces) <- "traces"
  .tracesTest(combi_traces)
  return(combi_traces)
}

#' Cluster peptides of gene to unique proteoforms across conditions
#' @param tracesList Object of class tracesList.
#' @param tracesListRaw Object of class tracesList. Default is NULL.
#' @param method Character string defining method for clustering
#' Default is "complete".
#' @param clusterH Numeric Cluster hight for cutree. Default is 0.5.
#' @param clusterN Integer Number of clusters. Default is NULL.
#' @param index Character string for method of cluster number estimation
#' Default is "silhouette".
#' @param plot logical,wether to print SibPepCorr density plot to R console.
#' Deafult is \code{FALSE}.
#' @param PDF logical, wether to print SibPepCorr density plot to a PDF file.
#' Deafult is \code{FALSE}.
#' @param name Character string with name of the plot, only used if
#' '\code{PDF=TRUE}.PDF file is saved under name.pdf.
#' Default is "hclustMultiCond".
#' @return Object of class tracesList with proteoform annotation
#' across conditions
#' @export
clustPepMultiCond <- function(tracesList, tracesListRaw=NULL,
                    method = c("complete","single"), clusterH = 0.5,
                    clusterN = NULL, index = "silhouette",
                    plot = FALSE, PDF=FALSE, name="hclustMultiCond", ...) {
  .tracesListTest(tracesList)
  method <- match.arg(method)
  if (!is.null(tracesListRaw)){
    .tracesListTest(tracesListRaw)
  } else {
    tracesListRaw <- tracesList
  }
  traces <- combineTracesMutiCond(tracesListRaw)
  clustered_traces <- clusterPeptides(traces,
                      method = method, clusterH = clusterH,
                      clusterN = clusterN, index = index,
                      plot = plot, PDF=PDF, name=name)
  clust_info <- subset(clustered_traces$trace_annotation,
                  select=c("id","protein_id","cluster","proteoform_id"))
  res <- lapply(tracesList, function(t){
   t$trace_annotation <- merge(t$trace_annotation,clust_info,
                          all.x=T,all.y=F,by=c("id","protein_id"))
   return(t)
  })
  class(res) <- "tracesList"
  .tracesListTest(res)
  return(res)
}

#' Annotate tracesList of multiple conditions with proteoforms that were
#' assigned for a combined traces object.
#' @param tracesList Object of class tracesList.
#' @param combinedTraces Object of class traces with annotated proteoforms.
#' @return Object of class tracesList witha nnotated proteoforms.
#' @export
annotateProteoformsAcrossConditions <- function(tracesList,combinedTraces){
  .tracesListTest(tracesList)
  .tracesTest(combinedTraces)
  if (!"proteoform_id" %in% names(combinedTraces$trace_annotation)){
    stop("CombinedTraces have not been annotated!")
  }
  res <- lapply(tracesList, function(t){
    cols <- names(t$trace_annotation)[which(names(t$trace_annotation) %in% names(combinedTraces$trace_annotation))]
    t$trace_annotation <- merge(t$trace_annotation,combinedTraces$trace_annotation,
                          all.x=T,all.y=F,by=cols,sort=F)
    return(t)
  })
  class(res) <- "tracesList"
  .tracesListTest(res)
  return(res)
}
