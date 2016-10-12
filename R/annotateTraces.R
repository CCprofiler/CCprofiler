#traces = pepTraces
#annotation_table <- "/Users/ibludau/Downloads/uniprot-proteome%3AUP000000625.txt"
#annotation_table <- read.table(annotation_table,sep = "\t",header = TRUE)
#names(annotation_table) <- c("protein_id","protein_mw")
#annotation_table$protein_mw <- as.numeric(gsub(",",".",annotation_table$protein_mw))
#write.table(annotation_table,"/Volumes/ibludau-1/SEC/Ecoli_data/Uniprot_mass_annotation.txt",sep="\t",col.names = TRUE,row.names = FALSE,quote = FALSE)

annotateTraces <- function(traces, annotation_table){
  annotation_table <- fread(annotation_table)
  new_annotation <- merge(traces$trace_annotation, annotation_table, by = "protein_id", all.x=TRUE)
  setcolorder(new_annotation, c("id", setdiff(names(new_annotation), "id")))
  traces$trace_annotation <- new_annotation
  traces
}