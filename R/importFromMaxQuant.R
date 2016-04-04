# Due to: http://stackoverflow.com/questions/24501245/data-table-throws-object-not-found-error
.datatable.aware=TRUE

#' Import peptide profiles from an MaxQuant DDA data analysis search result table.
#' @import data.table
#' @param file.name Quantitative MS data result file. Defaults to "peptides.txt"
#' @param annotation.table experimentalDesignTemplate table from MaxQuant, containing columns "Name", "Experiment" and "Fraction".
#'        Fraction should contain the fraction_number.
#' @param quanttype Type of quant to extract, defaults to "Intensity " and can be replaced by other prefixes as given in the MaxQuant result table header.
#'        E.g. "^Intensity H " or "^Intensity L " can be used to differentiate channels of SILAC-label data.
#' @return An object of class Traces (peptide level) containing \itemize{table.wide, table.long and id} that can be processed with the herein
#'         contained functions.
#'         
#' @export
#' 
importFromMaxQuant <- function(file.name = "peptides.txt", quanttype = "^Intensity ",
                                annotation.table = "experimentalDesignTemplate.txt") {
  # read data & annotation table & clean up
  ##################################################
  # data  <- data.table::fread(file.name, sep="\t", header=TRUE) #produces funny numbers for one of 88 columns (?)
  data <- as.data.table(read.csv(file.name, sep ="\t", header = TRUE, check.names = FALSE))
  message("reading results file ...")
  annotation <- data.table::fread(annotation.table)
  data.s <- data[grep("sp\\|", data$Proteins)] #retain only those peptides with a "sp|*" fasta entry in Proteins column
  
  # Select Intensity columns of proteotypic peptides 
  ##################################################
  Quantcolumns <- grep(quanttype, names(data.s))
  Labelcolumns <- which(names(data.s) %in% c("Sequence", "Proteins"))
  
  data.s.traces <- data.s[`Unique (Proteins)` == "yes", Quantcolumns, with = FALSE]
  data.s.labels <- data.s[`Unique (Proteins)` == "yes", Labelcolumns, with = FALSE]
  data.s.labels[, protein_id:=sapply(as.character(data.s.labels$Proteins), function(x){strsplit(x, split = "\\|")[[1]][2]})]
  data.s.labels[, protein_name:=sapply(as.character(data.s.labels$Proteins), function(x){strsplit(x, split = "\\|")[[1]][3]})]
  
  #Replace column headers by fraction number and order ascending
  nruns <- length(names(data.s.traces))
  column_headers_experiment <- gsub(quanttype, "", names(data.s.traces))
  column_headers_fraction <- sapply(column_headers_experiment, function(x){annotation[Experiment %in% x, "Fraction", with = FALSE]})
  column_headers_fraction <- sapply(column_headers_fraction, function(x){x[[1]]})
  names(data.s.traces) <- sprintf("%02d", column_headers_fraction)
  data.s.traces <- data.s.traces[,order(as.numeric(names(data.s.traces))), with = FALSE]
  
  # Assemble and output result "Traces" object
  #################################################
  traces.wide <- data.table(protein_id = data.s.labels$protein_id, peptide_id = data.s.labels$Sequence, data.s.traces)
  setorder(traces.wide, protein_id)
  traces.long <- melt(traces.wide,
                      id.vars = c("protein_id","peptide_id"),
                      value.name = "Intensity", 
                      variable.name = "fraction_number")
  labels <- data.s.labels[,c(2, 4, 3, 1), with = FALSE]
  names(labels)[4] <- "peptide_id"
  names(labels)[1] <- "fasta_id"
  setorder(labels, protein_id)
  result <- list(table.wide = traces.wide, table.long = traces.long, ids = labels)
  class(result) <- "Traces"
  return(result)
}
