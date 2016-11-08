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
#' @param also import nonproteotypic peptides corresponding to a single
#'      canonical uniprot id(default FALSE)
#' @return An object of class Traces containing 
#'     "traces", "traces_type", "traces_annotation" and "fraction_annotation" list entries
#'     that can be processed with the herein contained functions.
#' @export
importFromOpenSWATH <- function(data= 'OpenSwathData.tsv', 
                                annotation.table="annotation.txt",
                                remove_requantified=TRUE,
                                MS1Quant=FALSE, rm.decoy = FALSE,
                                collapseIsoforms = FALSE)
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

  
  # # convert ProteinName to uniprot ids
  # if (length(grep("\\|",data$ProteinName)) > 0) {
  #   message('converting ProteinName to uniprot ids ...')
  #   decoy_idx <- grep("^DECOY_",data$ProteinName)
  #   data$ProteinName <- gsub(".*\\|(.*?)\\|.*", "\\1", data$ProteinName)
  #   data$ProteinName[decoy_idx] <- paste0("DECOY_",data$ProteinName[decoy_idx])
  #   # the above does not wrk further downstream because "1/" is removed
  #   #data$ProteinName <- extractIdsFromFastaHeader(data$ProteinName)
  # }
  
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
  
  # Remove Decoys if specified
  
  if (rm.decoy == TRUE){
    message('removing decoys ...')
    n <- nrow(data.s)
    data.s <- data.s[!grep("^DECOY", data.s$ProteinName)]
    message('removed ', n - nrow(data.s), ' decoys')
  }
  
  # remove non-proteotypic peptides ,or peptides that match multiple canonical Protein Ids
  # (ifcollapseIsoforms == T)
  
  if (collapseIsoforms) {
    message("removing non-unique proteins only keeping peptides with single canonical Uniprot ID...")
    
    decoys_id <- grep("^DECOY", data.s$ProteinName)
    data.s$ProteinName[decoys_id] <- gsub("^DECOY_", "", data.s$ProteinName[decoys_id])
    isoform_no <- as.numeric(gsub("/.*", "", data.s$ProteinName))
    a <- strsplit(data.s$ProteinName, split = "\\|")
    uniprots <- lapply(1:length(a), function(i) a[[i]][seq(from = 2, to = isoform_no[i] * 2, by = 2)])
    # when dash is removed all ids should be equal
    canonical_typic <- sapply(uniprots, function(x) length(unique(gsub("-.*", "", x))) == 1)
    #Paste the full peptide name with all Isoforms onto the FullPeptideName
    data.s$FullPeptideName <- paste0(data.s$FullPeptideName, "-", data.s$ProteinName)
    data.s$ProteinName <- paste0(gsub("-.*", "", sapply(a, `[`, 2)))
    # data.s$ProteinName <- paste0(sapply(a,`[`,1),"|",
    # gsub("-.*","",sapply(a,`[`,2)),"|",gsub("/sp","",sapply(a,`[`,3)))
    data.s$ProteinName[decoys_id] <- paste0("DECOY_",data.s$ProteinName[decoys_id])
    data.s <- data.s[canonical_typic]
    message("Removed ", length(isoform_no) - nrow(data.s), " non unique peptides out of ",
            length(isoform_no))
    #data.s <- rbind(data.s, decoys)
    
  } else {
    
    #remove non-proteotypic discarding/keeping Decoys
    
    message('removing non-unique proteins only keeping proteotypic peptides ...')
    data.s <- data.s[c(grep("^1/", data.s$ProteinName), grep("^DECOY_1/", data.s$ProteinName))] 
    data.s$ProteinName <- gsub("1\\/","",data.s$ProteinName)
    
  }
  
  # Remove requantified peptides
  
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
