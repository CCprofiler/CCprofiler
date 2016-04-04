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
