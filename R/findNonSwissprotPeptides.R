#' Align peptides from a gene to the corresponding Swissprot sequence to extract
#' Peptides that would have been missed
#' @param pepTraces An object of type \code{traces.obj}.
#' @param proteinFasta Path to fasta (or \code{AAStringSet}) file with
#'proteinSequences used for alignment
#' @param idmapping Path to (or dataframe of) idMapping table obtainable from
#' ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/
#' @import Biostrings
#' @import data.table
#' @export

findNonSwissprotPeptides <- function(pepTraces, proteinFasta, id_table = NULL,
                                     mode = "ENSEMBL", download = FALSE, extr2 = T, alignsTo = F){
  # Read the data
  if(class(proteinFasta)[1] == "character"){
    sequences <- readAAStringSet(filepath = proteinFasta)
  }else if(class(proteinFasta)[1] == "AAStringSet"){
    sequences <- proteinFasta
  }else{
    stop("Input ProteinFasta must be of type character (Filepath to fasta file) or AAStringSet")
  }
  if(extr2){
    names <- strsplit(names(sequences), split = "\\|")
    names(sequences) <- sapply(names, "[", 2)
  }
  
  if(class(id_table)[1] == "character"){
    idmapping <- read.delim(id_table, 
                            sep = "\t", header = FALSE, stringsAsFactors = F)
    idmapping <- data.table(uniprot_id = idmapping[,"V1"], ensembl_id = idmapping[,"V19"])
    idmapping <- idmapping[, list(ensembl_id = unlist(strsplit(ensembl_id, ";"))), by=uniprot_id]
    
  }else if(any(class(id_table) %in% "data.table")){
    idmapping <- data.table(uniprot_id = id_table[,sp_id], ensembl_id = id_table[,from_id])
    idmapping <- idmapping[, list(ensembl_id = unlist(strsplit(ensembl_id, ";"))), by=uniprot_id]
  }else if(download){
    idmapping <- downloadIdMapping(pepTraces)
  }else{
    stop("Input id_table must be of type character(Path to table) or data.table")
  }
  setkey(idmapping, ensembl_id)
  
  # map the gene names
  
  alignsSwissprot <- apply(pepTraces$trace_annotation, 1, alignPeptide, sequences, idmapping, mode, alignsTo)
  pepTraces$trace_annotation$AlignsSwissprot <- alignsSwissprot
  return(pepTraces$trace_annotation)
}


alignPeptide <- function(Trace, sequences, idmapping, mode = "ENSEMBL", alignsTo){
  gene <- Trace[2]
  if(grepl("DECOY", gene)){
    return(FALSE)
  }else{
    pep <- gsub("\\(.*?\\)","",Trace[1])
    if(mode == "ENSEMBL"){
      up_ids <- idmapping[gene, uniprot_id]
    }else if (mode == "UNIPROT"){
      up_ids <- gene
    }
    # if(is.na(up_ids)) message(paste0("Warning: Gene ", gene, " had no corresponding uniprot ID in the mapping table"))
    seqs <- sequences[grep(up_ids, names(sequences))]
    if(length(seqs) == 0) message(paste0("Warning: Gene ", gene, 
                                         " had no corresponding Swissprot sequence"))
    matches <- grepl(pep, seqs)
    if(alignsTo){
      return(paste(names(seqs)[matches], collapse = ";"))
    }else{
      return(any(matches))
    }
  }
}

#' Download ENSG id to uniprotID mapping through Uniprot.ws
#' @param pepTraces An object of type \code{traces.obj}.
#' @import UniProt.ws
#' @import data.table
#' @export

downloadIdMapping <- function(pepTraces){
  gene_names <- unique(pepTraces$trace_annotation$protein_id)
  gene_names <- gene_names[!grepl("DECOY|sp", gene_names)]
  up <- UniProt.ws(taxId=9606)
  res <- select(up,
                keys = gene_names,
                columns = c("ENSEMBL", "UNIPROTKB"),
                keytype = "ENSEMBL")
  names(res) <- c("ensembl_id", "uniprot_id")
  res <- as.data.table(res)
  return(res)
  
}


