# Due to: http://stackoverflow.com/questions/24501245/data-table-throws-object-not-found-error
.datatable.aware=TRUE

#' Annotate traces object.
#' @description Add columns to the trace_annotation table from an external annotation table.
#' @import data.table
#' @param traces object of class traces
#' @param annotation_table tab-separated annotation table to be merged into trace_annotation, file or R data.table
#' @param traces_id_column name of the id column specifying the entries to be annotated
#' @param annotation_table_id_column name of the column in the reference table that contains the ids
#' defaults to "Entry" as expexted from uniprot.tab files
#' @param replace_whitespace whether spaces contained in the anntation table column names shall be replaced by "_"
#' defaults to TRUE
#' ,typically FullPeptideName or Uniprot Entry names as protein_id. Same as in traces$traces$id
#' @return traces object of class traces with extended annotation
#' @export

annotateTraces <- function(traces,
                           annotation_table = "uniprot.tab",
                           traces_id_column = "id",
                           annotation_table_id_column = "Entry",
                           replace_whitespace = TRUE){


  # use or if path to file read annotation table and add id column
  if (isTRUE(class(annotation_table) == "character")){
    annotation_table <- fread(annotation_table, header = TRUE)
  } else if (isTRUE(class(annotation_table) == "data.frame")){
    annotation_table <- as.data.table(annotation_table)
  }

  # test if the annotation column ids match those in traces
  if(sum(traces$trace_annotation[[traces_id_column]] %in% annotation_table[[annotation_table_id_column]]) < 1) {
    stop(paste("The ids in the annotation_table column", annotation_table_id_column,
    "do not match the IDs in the traces$trace_annotation$id column.")   )
  }

  if (replace_whitespace) {
   setnames(annotation_table,names(annotation_table),gsub(" ", "_", names(annotation_table)))
  }

  # If there's the Uniprot Mass column, extract numeric protein_mw in kDa
  if ("Mass" %in% names(annotation_table)){
    annotation_table[,protein_mw := as.numeric(gsub(",","",Mass))]
    annotation_table[,protein_mw := protein_mw/1000]
  }

  new_annotation <- merge(traces$trace_annotation,
                          annotation_table, by.x = traces_id_column, by.y = annotation_table_id_column , all.x=TRUE)
  traces$trace_annotation <- new_annotation

  # sort both trace_annotation and traces in case merge conducted an auto-sort
  setorder(traces$trace_annotation, id)
  setorder(traces$traces, id)
  return(traces)
}
