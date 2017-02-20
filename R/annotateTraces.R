#' Import peptide profiles from an OpenSWATH experiment.
#' @import data.table
#' @param traces.obj An object of type \code{traces.obj}.
#' @param annotation.table path to tab-separated .txt file for annotation of the traces
#' @return traces.obj An object of type \code{traces.obj}.
#' @export

#' annotate Traces object
#' @description Add columns to the trace_annotation table from an external annotation table
#' @param traces Object of class traces
#' @param annotation_table path to tab-separated annotation table to be merged into trace_annotation
#' @param id_column name/header of the column in the reference table that contains the ids
#' defaults to "Entry" as expexted from uniprot.tab files
#' @param replace_whitespace whether spaces contained in the anntation table column names shall be replaced by "_"
#' defaults to TRUE
#' ,typically FullPeptideName or Uniprot Entry names as protein_id. Same as in traces$traces$id
annotateTraces <- function(traces, annotation_table = "uniprot.tab", id_column = "Entry", replace_whitespace = TRUE){
  
  # use or if path to file read annotation table and add id column
  if (isTRUE(class(annotation_table) == "character")){
    annotation_table <- fread(annotation_table, header = TRUE)
  } else if (isTRUE(class(annotation_table) == "data.frame")){
    annotation_table <- as.data.table(annotation_table)
  }
  annotation_table$id = annotation_table[[id_column]]
  
  if (replace_whitespace) {
   names(annotation_table) <- gsub(" ", "_", names(annotation_table))
  }
  
  # If there's the Uniprot Mass column, extract numeric protein_mw in kDa
  if ("Mass" %in% names(annotation_table)){
    annotation_table[,protein_mw := as.numeric(gsub(",","",Mass))]
    annotation_table[,protein_mw := protein_mw/1000]
  }
  
  
  
  new_annotation <- merge(traces$trace_annotation,
                          annotation_table, by = "id", all.x=TRUE)
  traces$trace_annotation <- new_annotation
  
  # sort both trace_annotation and traces in case merge conducted an auto-sort
  setorder(traces$trace_annotation, id)
  setorder(traces$traces, id)
  return(traces)
}
