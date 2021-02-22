# Due to: http://stackoverflow.com/questions/24501245/data-table-throws-object-not-found-error
.datatable.aware=TRUE

#' Import peptide profiles from an OpenSWATH experiment.
#' @description This is a convenience function to directly import peptide profiles
#'  from an OpenSWATH experiment (after TRIC alignment). The peptide intensities are calculated by summing all charge
#'  states. Alternativley the MS1 signal can be used for quantification.
#' @import data.table
#' @param data Quantitative MS data in form of OpenSWATH result file or R data.table.
#' @param annotation_table file or data.table containing columns `filename` and
#'     `fraction_number` that map the file names (occuring in input table filename column)
#'     `to a SEC elution fraction.
#' @param rm_requantified Logical, whether requantified (noise) peak group quantities
#'     (as indicated by m_score = 2) should be removed, defaults to \code{TRUE}.
#' @param rm_decoys Logical, whether decoys should be removed, defaults to \code{FALSE}.
#' @param rm_nonProteotypic Logical, whether non-proteotypic peptides should be removed, defaults to \code{TRUE}.
#' @param proteogenomicsWF Logical, whether information on protein, isoform and gene level is encoded in ProteinName, defaults to \code{FALSE}.
#' @param MS1Quant Logical, whether MS1 quantification should be used, defaults to \code{FALSE}.
#' @param verbose Logical, whether to print progress message into console, defaults to \code{TRUE}.
#' @return An object of class Traces containing
#'     "traces", "traces_type", "traces_annotation" and "fraction_annotation" list entries
#'     that can be processed with the herein contained functions.
#' @export
#' @examples
#'   input_data <- exampleOpenSWATHinput
#'   annotation <- exampleFractionAnnotation
#'   traces <- importFromOpenSWATH(data = input_data,
#'                                 annotation_table = annotation,
#'                                 rm_requantified = TRUE)
#'   summary(traces)
#'
importFromOpenSWATH <- function(data,
                                annotation_table,
                                rm_requantified=TRUE,
                                rm_decoys = FALSE,
                                rm_nonProteotypic = TRUE,
                                MS1Quant=FALSE,
                                proteogenomicsWF=FALSE,
                                uniprot=TRUE,
                                verbose=TRUE){

  ## test arguments
  if (missing(data)){
        stop("Need to specify data in form of OpenSWATH result file or R data.table.")
  }
  if (missing(annotation_table)){
        stop("Need to specify annotation_table.")
  }

  ## read data & annotation table

  if (class(data)[1] == "character") {
    if (file.exists(data)) {
      message('reading results file ...')
      data  <- data.table::fread(data, header=TRUE)
    } else {
      stop("data file doesn't exist")
    }
  } else if (all(class(data) != c("data.table","data.frame"))) {
    stop("data input is neither file name or data.table")
  }

  if (class(annotation_table)[1] == "character") {
    if (file.exists(annotation_table)) {
      message('reading annotation table ...')
      annotation_table  <- data.table::fread(annotation_table)
    } else {
      stop("annotation_table file doesn't exist")
    }
  } else if (all(class(annotation_table) != c("data.table","data.frame"))) {
    stop("annotation_table input is neither file name or data.table")
  }

  ## discarding/keeping Decoys
  if (rm_decoys == TRUE){
    message('remove decoys ...')
    data <- data[grep("^DECOY", data$ProteinName,invert=TRUE)]
  }

  ## discarding/keeping non-proteotypic peptides
  if (rm_nonProteotypic ==TRUE) {
    message('remove non-proteotypic peptides ...')
    data <- data[grep("^1/|^DECOY_1/", data$ProteinName)]
    data$ProteinName <- gsub("1\\/","",data$ProteinName)
  }
  else {
    ## better to raise error than failing
    stop('No proteotypic peptides present in OSW file')
  }

  ## subset data to some important columns to save RAM
  if (MS1Quant == TRUE) {
  column_names <- c('transition_group_id', 'ProteinName','FullPeptideName',
                    'filename', 'Sequence', 'decoy', 'aggr_prec_Peak_Area',
                    'd_score', 'm_score')
  }else{
  column_names <- c('transition_group_id', 'ProteinName','FullPeptideName',
                      'filename', 'Sequence', 'decoy', 'd_score', 'm_score', 'Intensity')
  }
  data_s <- subset(data, select=column_names)

  ## Use aggregated precursor (MS1) area as Intensity column if selected
  if (MS1Quant == TRUE){
    setnames(data_s, 'aggr_prec_Peak_Area', 'Intensity')
  }

  rm(data)
  gc()

  ## remove requantified values if selected
  if (rm_requantified == TRUE) {
      data_s <- data_s[m_score < 2, ]
  }

  ## add fraction number column to main dataframe
  fraction_number <- integer(length=nrow(data_s))
  files <- annotation_table$filename
  data_filenames <- data_s$filename

  ## filter filenames to keep only present in annotation
  data_s <- data_s[grep(files, data_s$filename), ]

  for (i in seq_along(files)) {
      idxs <- grep(files[i], data_filenames)
      fraction_number[idxs] <- annotation_table$fraction_number[i]
      if (verbose) {
        message(paste("PROCESSED", i, "/", length(files), "filenames"))
      }
  }
  data_s <- cbind(data_s, fraction_number)

  ## Assemble and output result "traces" object
  traces_wide <-
      data.table(dcast(data_s, ProteinName + FullPeptideName ~ fraction_number,
                       value.var="Intensity",
                       ## averaging across charge states not summing
                       fun.aggregate=mean))

  traces_annotation <- data.table(traces_wide[,c("FullPeptideName", "ProteinName"), with = FALSE])
  setcolorder(traces_annotation, c("FullPeptideName", "ProteinName"))
  if (proteogenomicsWF==TRUE){
    message("Separating information on protein-, isoform- and gene-level within ProteinName...")
    setnames(traces_annotation,c("FullPeptideName"),c("id"))
    traces_annotation <- separateProteinNamesToGeneLevel(traces_annotation)

    traces_annotation[,decoy := 0]
    if(rm_decoys == FALSE){
      traces_annotation$decoy[grep("^DECOY", traces_annotation$ProteinName)] = 1
    }
    traces_annotation[,ProteinName:=NULL]
  } else {
    ## convert ProteinName to uniprot ids
    if (length(grep("\\|",traces_annotation$ProteinName)) > 0) {
      message('converting ProteinName to uniprot ids ...')
      decoy_idx <- grep("^DECOY_",traces_annotation$ProteinName)
      traces_annotation$ProteinName <- gsub(".*\\|(.*?)\\|.*", "\\1", traces_annotation$ProteinName)
      if(length(decoy_idx>0)) {
        traces_annotation$ProteinName[decoy_idx] <- paste0("DECOY_",traces_annotation$ProteinName[decoy_idx])
      }
    }
    setnames(traces_annotation,c("FullPeptideName","ProteinName"),c("id","protein_id"))
  }

  traces <- subset(traces_wide, select = -ProteinName)
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

  # sort both trace_annotation and traces for consistency (also done in annotateTraces function)
  setorder(result$trace_annotation, id)
  setorder(result$traces, id)
  .tracesTest(result)
  return(result)
}


separateProteinNamesToGeneLevel <- function(ann_table){
  split <- strsplit(ann_table$ProteinName, split = ";")
  newprot <- lapply(split, mapSinglePeptide)
  ann_table[,ensembl_protein_id:=unlist(newprot)[ c(TRUE,FALSE,FALSE) ]]
  ann_table[,isoform_id:=unlist(newprot)[ c(FALSE,TRUE,FALSE) ]]
  ann_table[,gene_id:=unlist(newprot)[ c(FALSE,FALSE,TRUE) ]]
  return(ann_table)
}

mapSinglePeptide <- function(proteins){
  idx_junc <- grep("^JUNC",proteins)
  if (length(idx_junc) > 0) {
    protein_map <- proteins[-idx_junc]
  } else {
    protein_map <- proteins
  }
  protein_map <- gsub("^cf\\|","",protein_map)
  isoform_map <- gsub("\\|.*$","",protein_map)
  gene_map <- unique(gsub("\\-.*$","",isoform_map))
  isoform_map <- unique(isoform_map)
  protein_map <- unique(gsub(".*\\|","",protein_map))
  isoform_map <- paste(c(length(isoform_map),isoform_map), collapse = "/")
  gene_map <- paste(c(length(gene_map),gene_map), collapse = "/")
  protein_map <- paste(c(length(protein_map),protein_map), collapse = "/")
  return(list(protein_map,isoform_map,gene_map))
}

#' Import peptide profiles from an OpenSWATH experiment.
#' @description This is a convenience function to directly import peptide profiles
#'  from an OpenSWATH experiment (after TRIC alignment). The peptide intensities are calculated by summing all charge
#'  states. Alternativley the MS1 signal can be used for quantification.
#' @import data.table
#' @param data Quantitative MS data in form of OpenSWATH result file or R data.table.
#' @param annotation_table file or data.table containing columns `filename` and
#'     `fraction_number` that map the file names (occuring in input table filename column)
#'     `to a SEC elution fraction.
#' @param rm_requantified Logical, whether requantified (noise) peak group quantities
#'     (as indicated by m_score = 2) should be removed, defaults to \code{TRUE}.
#' @param rm_decoys Logical, whether decoys should be removed, defaults to \code{FALSE}.
#' @param rm_nonProteotypic Logical, whether non-proteotypic peptides should be removed, defaults to \code{TRUE}.
#' @param proteogenomicsWF Logical, whether information on protein, isoform and gene level is encoded in ProteinName, defaults to \code{FALSE}.
#' @param MS1Quant Logical, whether MS1 quantification should be used, defaults to \code{FALSE}.
#' @param verbose Logical, whether to print progress message into console, defaults to \code{TRUE}.
#' @return An object of class Traces containing
#'     "traces", "traces_type", "traces_annotation" and "fraction_annotation" list entries
#'     that can be processed with the herein contained functions.
#' @export
#' @examples
#'   input_data <- exampleOpenSWATHinput
#'   annotation <- exampleFractionAnnotation
#'   traces <- importFromOpenSWATH(data = input_data,
#'                                 annotation_table = annotation,
#'                                 rm_requantified = TRUE)
#'   summary(traces)
#'
importMultipleCondiionsFromOpenSWATH <- function(data,
                                annotation_table,
                                rm_requantified=TRUE,
                                rm_decoys = FALSE,
                                rm_nonProteotypic = TRUE,
                                MS1Quant=FALSE,
                                proteogenomicsWF=FALSE,
                                verbose=TRUE){
  samples <- unique(annotation_table$sample)
  traces_list <- lapply(samples, function(s){
    ann_sub <- annotation_table[sample == s]
    dt_sub <- data[ann_sub[, filename]]
    print(paste("Importing: ", s))
    importFromOpenSWATH(data = dt_sub, annotation_table = ann_sub,
    rm_requantified=rm_requantified,
    rm_decoys = rm_decoys,
    rm_nonProteotypic = rm_nonProteotypic,
    MS1Quant=MS1Quant,
    proteogenomicsWF=proteogenomicsWF,
    verbose=verbose)
  })
  names(traces_list) <- samples
  class(traces_list) <- "tracesList"
  .tracesListTest(traces_list)
  traces_list
}
