#' Summarize a traces object
#' @export
summarize_proteogenomics <- function(traces){
  UseMethod("summarize_proteogenomics", traces)
}

#' Summarize a traces object
#' @export
summarize_proteogenomics.traces <- function(traces) {
  .tracesTest(traces)
  no_peptides <- length(unique(traces$trace_annotation$id))
  no_genes <- length(unique(traces$trace_annotation$gene_id))
  isoforms <- unique(traces$trace_annotation$isoform_id)
  isoforms <- gsub("^[0-9]+\\/","",isoforms)
  isoforms <- paste(isoforms,collapse="/")
  isoforms <- unique(unlist(strsplit(isoforms,split="/")))
  no_isoforms <- length(isoforms)
  no_isotypic <- length(unique(traces$trace_annotation$isoform_id[grep("^1/",traces$trace_annotation$isoform_id)]))
  res <- c(no_peptides, no_genes, no_isoforms, no_isotypic)
  names(res) <- c("No. of peptides", "No. of genes", "No. of isoforms", "No. of isotypic isoforms")

  gen_iso_map <- data.table(gene_id = traces$trace_annotation$gene_id,
                            isoform_id = gsub("^[0-9]+\\/","",traces$trace_annotation$isoform_id))
  gen_iso_map[,isoform_id := paste(isoform_id,collapse="/"), by=gene_id]
  gen_iso_map <- unique(gen_iso_map,by="gene_id")
  gen_iso_map[,isoform_id := paste(unique(unlist(strsplit(isoform_id,split="/"))),collapse="/"),by=gene_id]
  gen_iso_map[, n_isoforms := length(unlist(strsplit(isoform_id,split="/"))),by=gene_id]

  stat <- summary(gen_iso_map$n_isoforms)
  tab <- table(gen_iso_map$n_isoforms)

  summary <- list(ids=res,isoforms_per_gene_summary=stat,isoforms_per_gene_table=tab)
  summary
}

#' Summarize a tracesList object
#' @export

summarize_proteogenomics.tracesList <- function(traces) {
  .tracesListTest(traces)
  res <- lapply(names(traces),function(tr){
    trace <- traces[[tr]]
    cat(paste0("###########################\n## ",
           tr, "\n",
           "###########################\n\n"))
    print(summarize_proteogenomics.traces(trace))
  })

}


#' getProteogenomicsIds of traces
#' @export
getProteogenomicsIds <- function(traces){
  UseMethod("getProteogenomicsIds", traces)
}

#' getProteogenomicsIds of traces
#' @export
getProteogenomicsIds.traces <- function(traces){
  .tracesTest(traces)
  peptides <- unique(traces$trace_annotation$id)
  genes <- unique(traces$trace_annotation$gene_id)
  isoforms <- unique(traces$trace_annotation$isoform_id)
  isoforms <- gsub("^[0-9]+\\/","",isoforms)
  isoforms <- paste(isoforms,collapse="/")
  isoforms <- unique(unlist(strsplit(isoforms,split="/")))
  isotypic <- unique(traces$trace_annotation$isoform_id[grep("^1/",traces$trace_annotation$isoform_id)])
  list(peptides=peptides,genes=genes,isoforms=isoforms,isotypic=isotypic)
}

#' getProteogenomicsIds of tracesList object
#' @export
getProteogenomicsIds.tracesList <- function(tracesList){
  .tracesListTest(tracesList)
  res <- lapply(names(tracesList),function(tr){
    trace <- tracesList[[tr]]
    getProteogenomicsIds.traces(trace)
    })
  names(res) <- names(tracesList)
  res
}

#' Annotate peptides that are specific for a non major isoform of the gene of origin
#' @details Defines the most frequently detected isoform as the major isoform for that gene
#' and then detects peptides that do not match to that gene.
#' @export
annotateNonMajorIsoformPeptides <- function(traces){
  UseMethod("annotateNonMajorIsoformPeptides", traces)
}

#' @describeIn annotateNonMajorIsoformPeptides Annotate a single traces object
annotateNonMajorIsoformPeptides.traces <- function(traces){
  ann <- traces$trace_annotation
  isotable <- ann[, .(isoid= paste(ensembl_protein_id, collapse = "/")), by=gene_id]
  genes <- isotable$gene_id
  isoforms <- strsplit(isotable$isoid, split="/")
  names(isoforms) <- genes
  majorIso <- lapply(isoforms, function(gene){
    iso <- grep(gene, pattern="ENS", value=T)
    detections <- table(iso)
    major <- names(detections)[which.max(detections)]
    return(major)
  })
  ann[, majorIsoform := majorIso[gene_id]]
  ann[, nonMajorSpecific := !grepl(majorIsoform, ensembl_protein_id), by=id]
  traces$trace_annotation <- ann
  return(traces)

}

#' @describeIn annotateNonMajorIsoformPeptides Annotate a tracesList object
annotateNonMajorIsoformPeptides.tracesList <- function(tracesList){
  res <- lapply(names(tracesList), function(sample){
    message(paste("Annotating isoforms of sample", sample))
    traces <- tracesList[[sample]]
    annotateNonMajorIsoformPeptides.traces(traces = traces)
  })
  names(res) <- names(tracesList)
  class(res) <- "tracesList"
  return(res)
}
