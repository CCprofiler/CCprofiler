#TEST
#type = "UNIPROTKB"
#ids = c("P20248","BLA","Q9P1F3","Q12979","P09884","P28340","Q07864")
#type = "UNIPROTKB_name"
#ids=c("CCNA2_HUMAN","BLA","PAL4E_HUMAN","S23IP_HUMAN")
#type="Gene_names"
#ids=c("ABR","ABRACL","BLA","C6orf115")
#traces.obj <- traces.obj
#convertIDs(ids,type,traces.obj)

convertIDs <- function(ids=ids,type="UNIPROTKB",traces.obj=traces.obj){
  human_annotated_more <- uniprot_reviewed_9606_10082016_annotated
  if(type=="UNIPROTKB"){
    # which input names are not defined in our mapping table as UNIPROTKB
    not_defined <- ids[which(! ids %in% human_annotated_more$UNIPROTKB)]
    # remaining protein uniprot ids
    protein_ids=unique(human_annotated_more$UNIPROTKB[which(human_annotated_more$UNIPROTKB %in% ids)])
    # which input ids do not have an MS signal in the traces object
    no_ms_signal <- protein_ids[which(! protein_ids %in% traces.obj$traces$id)]
    # final list of uniprot accessions used for the analysis
    protein_ids <- protein_ids[which(protein_ids %in% traces.obj$traces$id)]
  } else if (type=="UNIPROTKB_name") {
    # which input names are not defined in our mapping table as UNIPROTKB_name
    not_defined <- ids[which(! ids %in% human_annotated_more$UNIPROTKB_name)]
    # remaining protein uniprot ids
    protein_ids=unique(human_annotated_more$UNIPROTKB[which(human_annotated_more$UNIPROTKB_name %in% ids)])
    # which input ids do not have an MS signal in the traces object
    no_ms_signal <- protein_ids[which(! protein_ids %in% traces.obj$traces$id)]
    # convert back to the input UNIPROTKB_name
    no_ms_signal <- unique(human_annotated_more$UNIPROTKB_name[which((human_annotated_more$UNIPROTKB %in% no_ms_signal) & (human_annotated_more$UNIPROTKB_name %in% ids))])
    # final list of uniprot accessions used for the analysis
    protein_ids <- protein_ids[which(protein_ids %in% traces.obj$traces$id)]
  } else if (type=="Gene_names") {
    # which input names are not defined in our mapping table as Gene_names
    not_defined <- ids[which(! ids %in% human_annotated_more$Gene_names)]
    # remaining protein uniprot ids
    protein_ids=unique(human_annotated_more$UNIPROTKB[which(human_annotated_more$Gene_names %in% ids)])
    # which input ids do not have an MS signal in the traces object
    no_ms_signal <- protein_ids[which(! protein_ids %in% traces.obj$traces$id)]
    # convert back to the input Gene_names
    no_ms_signal <- unique(human_annotated_more$Gene_names[which((human_annotated_more$UNIPROTKB %in% no_ms_signal) & (human_annotated_more$Gene_names %in% ids))])
    # final list of uniprot accessions used for the analysis
    protein_ids <- protein_ids[which(protein_ids %in% traces.obj$traces$id)]
  } else if (type=="ENSEMBL") {
    # which input names are not defined in our mapping table as ENSEMBL
    not_defined <- ids[which(! ids %in% human_annotated_more$ENSEMBL)]
    # remaining protein uniprot ids
    protein_ids=unique(human_annotated_more$UNIPROTKB[which(human_annotated_more$ENSEMBL %in% ids)])
    # which input ids do not have an MS signal in the traces object
    no_ms_signal <- protein_ids[which(! protein_ids %in% traces.obj$traces$id)]
    # convert back to the input ENSEMBL
    no_ms_signal <- unique(human_annotated_more$ENSEMBL[which((human_annotated_more$UNIPROTKB %in% no_ms_signal) & (human_annotated_more$ENSEMBL %in% ids))])
    # final list of uniprot accessions used for the analysis
    protein_ids <- protein_ids[which(protein_ids %in% traces.obj$traces$id)]
  } else if (type=="HGNC") {
    # which input names are not defined in our mapping table as HGNC
    not_defined <- ids[which(! ids %in% human_annotated_more$HGNC)]
    # remaining protein uniprot ids
    protein_ids=unique(human_annotated_more$UNIPROTKB[which(human_annotated_more$HGNC %in% ids)])
    # which input ids do not have an MS signal in the traces object
    no_ms_signal <- protein_ids[which(! protein_ids %in% traces.obj$traces$id)]
    # convert back to the input HGNC
    no_ms_signal <- unique(human_annotated_more$HGNC[which((human_annotated_more$UNIPROTKB %in% no_ms_signal) & (human_annotated_more$HGNC %in% ids))])
    # final list of uniprot accessions used for the analysis
    protein_ids <- protein_ids[which(protein_ids %in% traces.obj$traces$id)]
  }
  if(type=="UNIPROTKB"){
    mapping_table <- unique(subset(human_annotated_more,select=c("UNIPROTKB","name"),UNIPROTKB %in% protein_ids))
  } else {
    mapping_table <- unique(subset(human_annotated_more,select=c("UNIPROTKB","name",type),UNIPROTKB %in% protein_ids))
  }
  return(list(not_defined=not_defined, no_ms_signal=no_ms_signal, protein_ids=protein_ids,mapping_table=mapping_table))
}