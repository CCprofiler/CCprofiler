# Due to: http://stackoverflow.com/questions/24501245/data-table-throws-object-not-found-error
.datatable.aware=TRUE

#' Annotate traces object
#' @description Add custom annotation columns to the trace_annotation in a traces object from an external annotation table, e.g. UniProt.
#' @import data.table
#' @param traces Object of class traces.
#' @param annotation_table Tab-separated annotation table to be merged into trace_annotation, file or R data.table.
#' @param traces_id_column Character string with name of the id column specifying the entries to be annotated.
#' Defaults to "protein_id".
#' @param annotation_table_id_column Character string with name of the column in the reference table that contains the ids.
#' Defaults to "Entry" as expexted from uniprot.tab files.
#' @param annotation_table_mass_column Character string with name of the column in the reference table that contains the protein mass annotation in kDa. If \code{uniprot_mass_format = TRUE}, this column should be called "Mass" and is in Da (PCprofiler will do the conversion). Defaults to "Mass".
#' @param uniprot_mass_format Logical, whether mass annotation is in the format of uniprot ("Mass" column with molecular weight in Da and comma as thouthands separator).
#' @param replace_whitespace Logical, whether spaces contained in the anntation table column names shall be replaced by "_".
#' Defaults to \code{TRUE}.
#' @return Traces object of class traces with extended annotation.
#' @examples
#' # Load some example data:
#' inputTraces <- examplePeptideTracesUnannotated
#' inputAnnotation <- exampleTraceAnnotation
#' # Run the annotation:
#' annotatedTraces <- annotateTraces(
#' traces=inputTraces,
#' annotation_table=inputAnnotation,
#' traces_id_column = "protein_id",
#' annotation_table_id_column = "Entry")
# Inspect annotation result:
#' annotatedTraces$trace_annotation
#' @export


annotateTraces <- function(traces,
                           annotation_table,
                           traces_id_column = "protein_id",
                           annotation_table_id_column = "Entry",
                           annotation_table_mass_column = "Mass",
                           uniprot_mass_format = TRUE,
                           replace_whitespace = TRUE){


  ## annotation_table: If path to file is provided read annotation table, if data.table is provided use it directly.
  if (class(annotation_table)[1] == "character") {
    if (file.exists(annotation_table)) {
      message('reading annotation table ...')
      annotation_table  <- data.table::fread(annotation_table, header = TRUE)
    } else {
      stop("annotation_table file doesn't exist")
    }
  } else if (all(class(annotation_table) != c("data.table","data.frame"))) {
    stop("annotation_table input is neither file name or data.table")
  }

  ## test if the annotation column ids match those in traces
  if(sum(traces$trace_annotation[[traces_id_column]] %in% annotation_table[[annotation_table_id_column]]) < 1) {
    stop(paste("The ids in the annotation_table column", annotation_table_id_column,
    "do not match the IDs in the traces$trace_annotation$id column.")   )
  }

  if (replace_whitespace) {
   setnames(annotation_table,names(annotation_table),gsub(" ", "_", names(annotation_table)))
  }

  # If there's the Uniprot Mass column, extract numeric protein_mw in kDa
  if (uniprot_mass_format) {
    if ("Mass" %in% names(annotation_table)){
      annotation_table[,protein_mw := as.numeric(gsub(",","",Mass))]
      annotation_table[,protein_mw := protein_mw/1000]
    } else {
      message("Mass is not in uniprot format (mass in column called \"Mass\"). Please use direct uniprot annotation table or set uniprot_mass_format=FALSE and provide a annotation_table_mass_column which contains the protein mass in kDa. PCprofiler continues without mass annotation.")
    }
  } else {
    if (annotation_table_mass_column %in% names(annotation_table)){
      annotation_table[,protein_mw := get(annotation_table_mass_column)]
    } else {
      message("No protein mass annotation found, PCprofiler continues without mass annotation.")
    }
  }

  new_annotation <- merge(traces$trace_annotation,
                          annotation_table, by.x = traces_id_column, by.y = annotation_table_id_column , all.x=TRUE)
  traces$trace_annotation <- new_annotation

  # sort both trace_annotation and traces in case merge conducted an auto-sort
  setorder(traces$trace_annotation, id)
  setorder(traces$traces, id)
  return(traces)
}
