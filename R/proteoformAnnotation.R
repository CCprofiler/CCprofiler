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

filterByMaxCorr <- function(traces, cutoff = 0.85, plot = FALSE, PDF=FALSE, name="maxCorrHist") {
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


calculateMinCorr <- function(genePeptideList, plot = FALSE, PDF=FALSE, name="minCorrHist") {
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
  return(minCorrMatrices[])
}


fitDistr <- function(prot_names = NULL, mincorrtable, distr = "gamma", split_value = 0.2, plot = FALSE, PDF=FALSE, name="GammaFitted") {
  if (is.null(prot_names)) {
    mincorrs <- mincorrtable$mincorr
  } else {
    prot_idxs <- which(prot_names %in% mincorrtable$protein_id)
    mincorrs <- mincorrtable$mincorr[prot_idxs]
  }
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


calculateSplicePval <- function(prot_names= NULL, mincorrtable, distr_fitted,
                                adj.method = "fdr", plot = FALSE, PDF=FALSE, name="SplicePval") {
  if (is.null(prot_names)) {
    mincorrs <- mincorrtable$mincorr
  } else {
    prot_idxs <- which(prot_names %in% mincorrtable$protein_id)
    mincorrs <- mincorrtable$mincorr[prot_idxs]
  }

# here we use pgamma directly in the function for our gamma distr example, if we want to
# use other distr, we need a general way to represent all ditribution, e.g.,
# "p + distr abbr."
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
  mincorrtable[, adj_pval := pval_adj]
  return(mincorrtable[])
}

clusterPeptides <- function(genePeptideList, method = "single", clusterH = 0.5, clusterN = NULL, index = "silhouette", plot = FALSE, PDF=FALSE, name="hclust", ...) {
  if (PDF) {
    pdf(paste0(name,".pdf"))
  }
  clustMatricesMethod <- lapply(genePeptideList, function(gene){
    genecorr <- cor(gene)
    if (is.null(clusterH)) {
      if (is.null(clusterN)) {
        nbClust_res <- NbClust(data=genecorr, diss=as.dist(1-genecorr), distance = NULL, min.nc=2, max.nc=min(4,nrow(genecorr)-1)[1],method = method, index = index)
        clusterN <- nbClust_res$Best.nc[1]
      }
    }
    cl <- hclust(as.dist(1-genecorr), method)
    if (plot == TRUE) {
      #plot(cl, labels = FALSE)
      plot(as.dendrogram(cl), ylim = c(0,1))
    }
    groups <- cutree(cl, k = clusterN, h = clusterH)
    groupsDT <- data.table(id=names(groups),cluster=groups)
    return(groupsDT)
  })
  if (PDF) {
    dev.off()
  }
  return(clustMatricesMethod)
}
