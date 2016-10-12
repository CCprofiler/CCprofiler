#' extractIdsFromFastaHeader
#' @description sp|P11021|GRP78_HUMAN --> P11021
#' @param fasta_ids character vector of fasta protein ids
#' @param id_type "accession" or "name"
#' @param keep.prefix whether IDs shall stay merged to the prefix 1/ or DECOY_1/
#' @export
extractIdsFromFastaHeader <- function(fasta_ids, id_type = "accession", keep.prefix = TRUE) {
  split <- strsplit(fasta_ids, split = "\\|")
  prefix <- up_accession <- sapply(split, function(x){x[[1]][1]})
  prefix <- gsub("sp", "", prefix)
  up_accession <- sapply(split, function(x){x[[2]][1]})
  up_name <- sapply(split, function(x){x[[3]][1]})
  if (keep.prefix){
    up_accession <- paste0(prefix, up_accession)
    up_name <- paste0(prefix, up_name)
  }
  if (id_type == "accession"){
    return(up_accession)
  }
  if (id_type == "name"){
    return(up_name)
  }
}

