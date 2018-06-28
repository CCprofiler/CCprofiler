
#' Find the leading Isoform for each peptide
#' @param traces Object of class tracesList
#' @param isoform_col Character, name of the column in the =trace_annotation= table
#' containing the annotation of all found isoforms for that peptide. Note:
#' This function assumes isoform annotation in the format:
#' '3/ISOFORM1/ISOFORM2/ISOFORM3'
#' @param ouput_col character, name of the column that contains the leading isoform
#' @details This function goes through every gene and applies a maximum parsimony
#' approach to find the leading Isoform for each peptide. For ever gene
#' (defined either by the 'gene_id' or the 'protein_id' column)
#' and counts the number of occurences for every isoform. Peptides are then annotated
#' with the fewest possible number of isoforms that are sufficient to explain all peptides.
#' In case of ties one isoform is randomly selected.
#' @return A traces object with an additional column in =trace_annotation= containing
#' the leading isoform.
#' @export

annotateLeadingIsoform <- function(traces, isoform_col="isoform_id", output_col="LeadingIsoform"){
  UseMethod("annotateLeadingIsoform", traces)
  }

#' @describeIn annotateLeadingIsoform Annotate the Leading Isoform in a single traces object
annotateLeadingIsoform.tracesList <- function(tracesList, isoform_col="isoform_id",
                                              output_col="LeadingIsoform"){

  .tracesListTest(tracesList)
  res <- lapply(tracesList, annotateLeadingIsoform.traces,
                isoform_col=isoform_col,
                output_col=output_col)
  class(res) <- "tracesList"
  .tracesListTest(res)
  return(res)
}

#' @describeIn annotateLeadingIsoform Annotate the Leading Isoform in a single traces object
annotateLeadingIsoform.traces <- function(traces, isoform_col="isoform_id", output_col="LeadingIsoform"){
  .tracesTest(traces, type = "peptide")
  ann <- copy(traces$trace_annotation)
  if(!("gene_id" %in% colnames(ann))){
    message("No gene_id column found for grouping genes. Trying to use protein_id instead")
    if(!("protein_id" %in% colnames(ann))){
      stop("Column protein_id not found for grouping genes.
 Either 'gene_id' or 'protein_id' column must be present! Aborting.")
    }else{
      groups <- "protein_id"
    }
  }else{
    groups <- "gene_id"
  }

  if(!(isoform_col %in% colnames(ann))){
    stop(paste0("Could not find column: ", isoform_col))
  }
  setnames(ann, isoform_col, "IsoformContainingColumn")
  ann[[groups]] <- gsub("^.*?\\/", "", ann[[groups]])
  suppressWarnings(ann[,leadingIsoform := {
    iso <- strsplit(gsub("^.*?\\/", "", IsoformContainingColumn), split="/")
    names(iso) <- id
    peps <- do.call(c,lapply(id, function(id) rep(id,length(iso[[id]]))))
    isos <- do.call(c, iso)
    presenceMat <- table(peps, isos)
    findLeadingIso(presenceMat)
    ## iso <- lapply(iso, "[", )
  } , by= get(groups)])

  setnames(ann, "IsoformContainingColumn", isoform_col)
  setnames(ann, "leadingIsoform", output_col)
  traces$trace_annotation <- ann[order(traces$traces$id)]
  .tracesTest(traces)
  return(traces)
}

#' Convert a presence matrix of isoform annotations to a vector of leading Isoforms
#' @param presenceMat numeric matrix with peptides as rows and isoforms as columns,
#' where a one indicates that that peptide is part of an isoform and a 0 otherwhise.
#' @details Applies a maximum parsimony approach, where the least number of isoforms
#' are used to explain a peptide.
#' @return character vector with length = nrow(presenceMat), with the leading
#' Isoform of every peptide.

findLeadingIso <- function(presenceMat){
  leadingIso <- character(nrow(presenceMat))
  isos <- colnames(presenceMat)
  presenceMat <- presenceMat == 1
  s <- colSums(presenceMat)
  for(i in 1:length(s)){
    l <- which.max(s) # This choses a random (first) isoform when a tie is reached
    leadingIso[presenceMat[,l]] <- isos[l]
    if(all(leadingIso != "")) break
    presenceMat[presenceMat[,l],] <- F
    presenceMat <- presenceMat[,-l, drop=F]
    s <- s[-l]
    isos <- isos[-l]
  }
  return(leadingIso)
}
