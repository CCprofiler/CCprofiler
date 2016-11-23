# Due to: http://stackoverflow.com/questions/24501245/data-table-throws-object-not-found-error
.datatable.aware=TRUE

#' Import peptide profiles from an OpenSWATH experiment.
#' @import data.table
#' @param data Quantitative MS data in form of OpenSWATH result file or R data.table.
#' @param annotation.table path to tab-separated .txt file containing columns `filename` and
#'     `fraction_number` that map the file names (occuring in input table filename column)
#'     `to a SEC elution fraction.
#' @param remove_requantified Whether requantified (noise) peak group quantities
#'     (as indicated by m_score = 2) should be removed, defaults to TRUE.
#' @return An object of class Traces containing 
#'     "traces", "traces_type", "traces_annotation" and "fraction_annotation" list entries
#'     that can be processed with the herein contained functions.
#' @export
importFromOpenSWATH <- function(data= 'OpenSwathData.tsv', 
                                annotation.table="annotation.txt",
                                remove_requantified=TRUE,
                                MS1Quant=FALSE, rm.decoy = FALSE)
  {
  
  # read data & annotation table
  ##################################################
  
  if (class(data)[1] == "character") {
    message('reading results file ...')
    data  <- data.table::fread(data, sep="\t", header=TRUE)
  } else if (all(class(data) != c("data.table","data.frame"))) {
    stop("data input is neither file name or data.table")
  }
  
  annotation <- data.table::fread(annotation.table)
  
  #remove non-proteotypic discarding/keeping Decoys
  message('removing non-unique proteins only keeping proteotypic peptides ...')
  if (rm.decoy == TRUE){
    message('remove decoys ...')
    data <- data[grep("^1/", data$ProteinName)] 
  } else {
    data <- data[c(grep("^1/", data$ProteinName), grep("^DECOY_1/", data$ProteinName))] 
  }
  data$ProteinName <- gsub("1\\/","",data$ProteinName) 
  
  
  # convert ProteinName to uniprot ids
  if (length(grep("\\|",data$ProteinName)) > 0) {
    message('converting ProteinName to uniprot ids ...')
    decoy_idx <- grep("^DECOY_",data$ProteinName)
    data$ProteinName <- gsub(".*\\|(.*?)\\|.*", "\\1", data$ProteinName)
    data$ProteinName[decoy_idx] <- paste0("DECOY_",data$ProteinName[decoy_idx])
    # the above does not wrk further downstream because "1/" is removed
    #data$ProteinName <- extractIdsFromFastaHeader(data$ProteinName)
  }
  
  # subset data to some important columns to save RAM
  if (MS1Quant == TRUE) {
  column.names <- c('transition_group_id', 'ProteinName','FullPeptideName',
                    'filename', 'Sequence', 'decoy', 'aggr_prec_Peak_Area',
                    'd_score', 'm_score')
  }else{
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
      data.s <- data.s[m_score < 2, ]
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
  setcolorder(traces_annotation, c("FullPeptideName", "ProteinName"))
  setnames(traces_annotation,c("FullPeptideName", "ProteinName"),c("id","protein_id"))
  
  traces <- subset(traces.wide, select = -ProteinName)
  traces[,id:=FullPeptideName]
  traces[,FullPeptideName:=NULL]
  
  nfractions <- ncol(traces)-1
  fractions <- as.numeric(c(1:nfractions))
  fraction_annotation <- data.table(id=fractions)
  
  traces_type = "peptide"
  
  result <- list("traces" = traces,
                 "trace_type" = traces_type,
                 "trace_annotation" = traces_annotation,
                 "fraction_annotation" = fraction_annotation)
  class(result) <- "traces"
  names(result) <- c("traces", "trace_type", "trace_annotation", "fraction_annotation")

  return(result)
}
