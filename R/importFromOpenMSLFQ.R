# Due to: http://stackoverflow.com/questions/24501245/data-table-throws-object-not-found-error
.datatable.aware=TRUE

#' Import peptide profiles from an OpenMS MS1-LFQ DDA data analysis result table.
#' @import data.table
#' @param file.name Quantitative MS data result file. Defaults to "OpenMS_LFQ_peptides.csv"
#' @param annotation.table Annotation table containing columns "Name" (-> raw file name), and "Fraction".
#'        Fraction should contain the fraction_number.
#' @return An object of class Traces (peptide level) containing \itemize{table.wide, table.long and id} that can be processed with the herein
#'         contained functions.
#'         
#' @export
#' 
importFromOpenMSLFQ <- function(file.name = "OpenMS_LFQ_peptides.csv", annotation.table = "OpenMS_LFQ_annotation.txt") {
  
  # read data & annotation table & clean up
  ##################################################
  # data  <- data.table::fread(file.name, sep="\t", header=TRUE) #produces funny numbers for one of 88 columns (?)
  data <- fread(file.name, header = TRUE, integer64 = "double")
  annotation <- data.table::fread(annotation.table)
  
  # Select Intensity columns of proteotypic peptides 
  ##################################################
  data.s <- data[n_proteins == 1] #remove ambiguos peps
  Quantcolumns <- grep("~", names(data.s))
  Labelcolumns <- which(names(data.s) %in% c("peptide", "protein"))
  
  data.s.traces <- data.s[, Quantcolumns, with = FALSE]
  data.s.labels <- data.s[, Labelcolumns, with = FALSE]
  
  #Replace column headers by fraction number and order ascending
  headers <- names(data.s.traces)
  nruns <- length(names(data.s.traces))
  filenames <- sapply(headers, function(x){strsplit(x, split = "~")[[1]][1]})
  annotation$Name <- toupper(annotation$Name)
  fractions <- sapply(filenames, function(x){annotation[Name %in% x, "Fraction", with = FALSE]})
  fractions <- sapply(fractions, function(x){x[[1]]})
  names(data.s.traces) <- sprintf("%02d", fractions)
  data.s.traces <- data.s.traces[,order(as.numeric(names(data.s.traces))), with = FALSE]
  
  # Assemble and output result "Traces" object
  #################################################
  traces.wide <- data.table(protein_id = data.s.labels$protein, peptide_id = data.s.labels$peptide, data.s.traces)
  traces.wide <- traces.wide[grep("DECOY", traces.wide$protein_id, invert = TRUE)]
  #Revert to clean Uniprot ids, save fasta and name entry information separately
  fasta_id <- traces.wide$protein_id
  protein_id <- sapply(as.character(traces.wide$protein_id), function(x){strsplit(x, split = "\\|")[[1]][2]})
  protein_name <- sapply(as.character(traces.wide$protein_id), function(x){strsplit(x, split = "\\|")[[1]][3]})
  traces.wide$protein_id <- protein_id
  setorder(traces.wide, protein_id)
  
  traces.long <- melt(traces.wide,
                      id.vars = c("protein_id","peptide_id"),
                      value.name = "Intensity", 
                      variable.name = "fraction_number")
  labels <- data.table(fasta_id = fasta_id, protein_id = protein_id, protein_name = protein_name, peptide_id = traces.wide$peptide_id)
  setorder(labels, protein_id)
  
  result <- list(table.wide = traces.wide, table.long = traces.long, ids = labels)
  class(result) <- "Traces"
  return(result)
}