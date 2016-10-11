# For some pretty weird reason
.datatable.aware=TRUE


EnrichmentAnalyzer <- function(Gene_list, Ontology_table , ID_Term_mapping , mode = 'file') { #mode can be either file or object
  
  if (mode == 'file') {
    
    genes <- as.data.frame(data.table::fread(Gene_list, sep="\t", header=TRUE))
    amigo <- as.data.frame(data.table::fread(Ontology_table, sep="\t", header=TRUE))
    map <- as.data.frame(data.table::fread(ID_Term_mapping, sep="\t", header=TRUE))
    
  } else if (mode == 'object'){
    
    genes <- as.data.frame(Gene_list)
    amigo <- as.data.frame(Ontology_table)
    map <- as.data.frame(ID_Term_mapping)
    
  }
  
  GO <- merge(amigo,map,by = 'GO_ID',all.x = TRUE)
  GO_Term <- unique(as.vector(t(GO$GO_Term)))
  genes <- unique(as.vector(t(genes)))
  
  n_GO_genes <- length(unique(as.vector(GO$Uniprot_ID)))
  
  df_output <- as.data.frame(matrix(nrow = length(GO_Term), ncol = 6))
  colnames(df_output) <- c('GO_Term', 'Size', 'Status', 'Representative_value', 'Pvalue', 'Adj_Pvalue')
  
  for (i in 1:length(GO_Term)) {
    handle.term <- GO_Term[i]
    handle <- subset(GO, GO$GO_Term == handle.term)
    handle.status <- handle$Status[1]
    handle.genes <- unique(as.vector(handle$Uniprot_ID))
    mutual <- intersect(genes , handle.genes)
    
    pvalue <- phyper(length(mutual), length(handle.genes), n_GO_genes - length(handle.genes), length(genes), lower.tail = FALSE)
    
    rep_value <-  length(mutual) / (length(genes) * length(handle.genes) / n_GO_genes)
    
    df_output$GO_Term[i] <- GO_Term[i]
    df_output$Size[i] <- length(handle.genes)
    df_output$Status[i] <- handle.status
    df_output$Representative_value[i] <- rep_value
    df_output$Pvalue[i] <- pvalue
    
    rm (handle,handle.status,handle.genes,mutual,rep_value,pvalue)
    
  }
  
  
  df_output$Adj_Pvalue <- p.adjust(df_output$Pvalue, method = 'BH')
  
  return(df_output)
   
}