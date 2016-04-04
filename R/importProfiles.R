# Due to: http://stackoverflow.com/questions/24501245/data-table-throws-object-not-found-error
.datatable.aware=TRUE

#' Import peptide profiles from an OpenSWATH experiment.
#' @import data.table
#' @param file.name Quantitative MS data result file.
#' @param annotation.table A data frame with two columns `run_id` and
#'        `fraction_number` that map the run identifier as contained in `file.name`
#'        to a SEC elution fraction.
#' @param remove_requantified Whether requantified (noise) peak group quantities should be removed, defaults to TRUE.
#'        
#' @return An object of class Traces (peptide level) containing \itemize{table.wide, table.long and id} that can be processed with the herein
#'         contained functions.
#'         
#' @export
importFromOpenSWATH <- function(file.name = "OpenSwathData.tsv", 
                                annotation.table = "annotation.txt",
                                remove_requantified = TRUE) {
  # read data & annotation table
  ##################################################
  data  <- data.table::fread(file.name, sep="\t", header=TRUE)
  message('reading results file ...')
  annotation <- data.table::fread(annotation.table)
  column.names <- c('transition_group_id', 'ProteinName','FullPeptideName',
                    'filename', 'Sequence', 'decoy', 'aggr_prec_Peak_Area',
                    'd_score', 'm_score', 'Intensity')
  data.s <- subset(data, select = column.names)
  if (remove_requantified == TRUE){
    data.s <- data.s[m_score <= 1,]
  }
  
  # add fraction number column to main dataframe
  ##################################################
  fraction_number <- integer(length = nrow(data.s))
  files <- annotation$filename
  data.filenames <- data.s$filename
  
  if (length(files) != length(unique(data.filenames))) {
    stop("Number of file names in annotation does not match data")
    }
  
  for (i in seq_along(files)) {
    idxs <- grep(files[i], data.filenames)
    fraction_number[idxs] <- annotation$fraction_number[i]
    message(paste("PROCESSED", i, "/", length(files), "filenames"))
  }
  data.s <-cbind(data.s, fraction_number)
  
  # Assemble and output result "Traces" object
  #################################################
  traces.wide <- data.table(dcast(data.s, ProteinName + FullPeptideName ~ fraction_number, value.var = "Intensity", fun.aggregate = sum))
  names(traces.wide)[1:2] <- c("protein_id", "peptide_id")
  traces.long <- melt(traces.wide,
                      id.vars = c("protein_id","peptide_id"),
                    value.name = "Intensity", 
                    variable.name = "fraction_number")
  result <- list(table.wide = traces.wide, table.long = traces.long, ids = traces.wide[,1:2, with = FALSE])
  class(result) <- "Traces"
  return(result)
}


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