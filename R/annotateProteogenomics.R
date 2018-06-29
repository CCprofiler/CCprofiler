
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

#' @describeIn annotateLeadingIsoform Annotate the Leading Isoform in a tracesList object
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

#' Align peptides to reference proteins to get the relative position
#' @param traces Object of class traces with annotated leading isoform (see =annotateLeadingIsoform=)
#' @param mappingTable data.table containing all sequences that were in the proteomics search space.
#' Must have the following columns: 'sequence', 'header'
#' 'header' must be of the format: 'sp|IsoformId|etc'. The function will extract the isoform id
#' and match it to the leading isoform.
#' @param multimatch character, Skip or choose the first sequence if multiple id's match.
#' @param verbose boolean, whether to be verbose.
#' @return Object of class traces with annotated relative positions of each peptide.
#' Two columns are added to =trace_annotation=:
#' 'PeptidePositionStart' the position of the first amino acid in the protein.
#' 'PeptidePositionEnd' the position of the last amino acid in the protein.
#' @export

annotateRelativePepPos <- function(traces, mappingTable, multimatch=c("first", "skip"), verbose=F){
  UseMethod("annotateRelativePepPos", traces)
  }

#' @describeIn annotateRelativePepPos Annotate the relative Position of all peptides in a tracesList object.

annotateRelativePepPos.tracesList <- function(tracesList, mappingTable, multimatch=c("first", "skip"), verbose=F){
  .tracesListTest(tracesList)
  res <- lapply(names(tracesList),function(x){
    message(paste0("Annotating Peptide positions in ", x))
    annotateRelativePepPos.traces(tracesList[[x]],
                                  mappingTable=mappingTable,
                                  multimatch=multimatch,
                                  verbose=verbose)
  })
  names(res) <- names(tracesList)
  class(res) <- "tracesList"
  .tracesListTest(res)
  return(res)
}

#' @describeIn annotateRelativePepPos Annotate the relative Position of all peptides in a single traces object.

annotateRelativePepPos.traces <- function(traces, mappingTable, multimatch=c("first", "skip"), verbose=F){
  multimatch <- match.arg(multimatch)
  ann <- traces$trace_annotation
  ## Remove PTMs
  if(!("Sequence" %in% names(ann))){
    ann$Sequence <- gsub("\\(.*?\\)", "", ann$id)
  }
  mappingTable$IsoformId <- gsub("\\|.*", "", gsub(">.*?\\|","", mappingTable$header))
  ## Find protein and align
  matches <- apply(ann, 1, function(pep){
    seq <- mappingTable[IsoformId == pep["LeadingIsoform"]]$sequence
    if(length(seq) == 1){
      return(gregexpr(pep["Sequence"], seq))
    }else if(length(seq) > 1){
      if(multimatch == "skip"){
        if(verbose){
          message(paste0("Warning: Multiple sequences found for Protein: ",
                         pep["LeadingIsoform"], ". Skipping"))
        }
        m <- -1
        attr(m, "match.length") <- 0
        return(list(m))
      }else{
        if(verbose){
          message(paste0("Warning: Multiple sequences found for Protein: ",
                         pep["LeadingIsoform"], ". Choosing First Match."))
        }
        seq <- seq[1]
        return(gregexpr(pep["Sequence"], seq))
      }
    }else{
      if(verbose){
        message(paste0("Warning: No sequence found for Protein: ",
                       pep["LeadingIsoform"], ". Skipping"))
      }
      m <- -1
      attr(m, "match.length") <- 0
      return(list(m))
    }
  })

  ann$PeptidePositionStart <- sapply(matches, function(m) m[[1]][1])
  ann$PeptidePositionEnd <- ann$PeptidePositionStart +
    sapply(matches, function(m) attr(m[[1]], "match.length")[1]) - 1
  traces$trace_annotation <- ann
  return(traces)

}
