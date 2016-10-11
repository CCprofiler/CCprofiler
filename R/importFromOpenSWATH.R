# Due to: http://stackoverflow.com/questions/24501245/data-table-throws-object-not-found-error
.datatable.aware=TRUE

#' Import peptide profiles from an OpenSWATH experiment.
#' @import data.table
#' @param file.name Quantitative MS data result file.
#' @param annotation.table path to tab-separated .txt file containing columns `filename` and
#'     `fraction_number` that map the file names (occuring in input table filename column)
#'     `to a SEC elution fraction.
#' @param remove_requantified Whether requantified (noise) peak group quantities
#'     (as indicated by m_score = 2) should be removed, defaults to TRUE.
#' @return An object of class Traces containing 
#'     "traces", "traces_type", "traces_annotation" and "fraction_annotation" list entries
#'     that can be processed with the herein contained functions.
#' @export
importFromOpenSWATH <- function(file_name= 'OpenSwathData.tsv', mode = 'file', #either 'file' or 'object' 
                                annotation.table="annotation.txt",
                                remove_requantified=TRUE,
                                MS1Quant=FALSE, rm.decoy = FALSE)
  {
  
  # read data & annotation table
  ##################################################
  
  if (mode == 'file') {
    
    message('reading results file ...')
    data  <- data.table::fread(file_name, sep="\t", header=TRUE)
    
  } else if (mode == 'object'){
      data <- as.data.table(file_name)
    }
  
  
  
  #remove non-proteotypic discarding/keeping Decoys
  if (rm.decoy == TRUE){
    data <- data[grep("^1/", data$ProteinName)] 
  } else {
    data <- data[c(grep("^1/", data$ProteinName), grep("^DECOY_1/", data$ProteinName))] 
  }
  # data$ProteinName <- gsub("^1/","", data$ProteinName) #remove 2/ etc. tags
  annotation <- data.table::fread(annotation.table)
  
  # subset data to some important columns to save RAM
  if (MS1Quant == TRUE) {
  column.names <- c('transition_group_id', 'ProteinName','FullPeptideName',
                    'filename', 'Sequence', 'decoy', 'aggr_prec_Peak_Area',
                    'd_score', 'm_score')
  }
  else{
  column.names <- c('transition_group_id', 'ProteinName','FullPeptideName',
                      'filename', 'Sequence', 'decoy', 'd_score', 'm_score', 'Intensity')
  }
  
  data.s <- subset(data, select=column.names)
  
  # Use aggregated precursor (MS1) area as Intensity column if selected
  if (MS1Quant == TRUE){
    setnames(data.s, 'aggr_prec_Peak_Area', 'Intensity')
  }
  
  rm(data)
  gc()
  
  if (remove_requantified == TRUE) {
      data.s <- data.s[m_score <= 1, ]
  }
  
  # add fraction number column to main dataframe
  ##################################################
  fraction_number <- integer(length=nrow(data.s))
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
  data.s <- cbind(data.s, fraction_number)
  
  # Assemble and output result "Traces" object
  #################################################
  traces.wide <-
      data.table(dcast(data.s, ProteinName + FullPeptideName ~ fraction_number,
                       value.var="Intensity",
                       fun.aggregate=sum))
  
  traces_annotation <- data.table(traces.wide[,c("FullPeptideName", "ProteinName"), with = FALSE])
  traces_annotation[, id:=FullPeptideName]
  
  traces <- subset(traces.wide, select = -ProteinName)
  setnames(traces, "FullPeptideName", "id")
  
  traces_type = "peptide"
  
  nfractions <- length(names(traces))
  fractions <- as.numeric(names(traces)[2:nfractions])
  fraction_annotation <- data.table(fractions)
  
  traces <- as.data.frame(traces)
  traces <- traces[,c(2:nfractions,1)]
  traces <- as.data.table(traces)
  
  result <- list("traces" = traces,
                 "traces_type" = traces_type,
                 "traces_annotation" = traces_annotation,
                 "fraction_annotation" = fraction_annotation)
  class(result) <- "Traces"
  names(result) <- c("traces", "traces_type", "traces_annotation", "fraction_annotation")

  return(result)
}
