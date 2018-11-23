# Due to: http://stackoverflow.com/questions/24501245/data-table-throws-object-not-found-error
.datatable.aware=TRUE

#' Protein Quantification
#' @description Calculate protein quantities basen on the topN peptide intensities.
#' @import data.table
#' @param traces An object of type traces, trace_type must be peptide.
#' @param topN Numeric integer, specifying the number of peptides to sum for protein quantification. Default is 2.
#' @param keep_less Logical, specifying whether proteins with less than topN peptides
#'  should be kept in the data (This may result in some protein intensities being calculated
#'  as the sum of fewer peptide intensities than others. Only use withh caution.), default is \code{FALSE}.
#' @param rm_decoys Logical, specifying whether decoys should be kept.
#'  The decoys have only limited use on the protein level, default is \code{TRUE}.
#' @param use_sibPepCorr Logical, specify whether the sibling peptide correlation should be considered for selecting
#' the top N peptides, default is \code{FALSE}.
#' @param use_repPepCorr Logical, specify whether the replicate peptide correlation should be considered for selecting
#' the top N peptides, default is \code{FALSE}. This argument is only used for tracesList objects.
#' @param full_intersect_only Logical, specify if only peptides detected across all traces in tracesList object should
#' should be used for quantification. Be aware that this might significantly reduce the number of quantifyable proteins.
#' This argument is only used for tracesList objects.
#' @param verbose Logical if warning messages should be printed. Default is TRUE.
#' @return An object of type traces, trace_type is protein.
#' @export
#' @examples
#' ## Load example data
#' pepTraces <- examplePeptideTracesFiltered
#'
#' ## Sum the intensities of the top 2 peptides to get protein intensities
#' protTraces <- proteinQuantification(pepTraces,
#'                                      topN = 2)
#' ## Check the result
#' summary(pepTraces)
#' summary(protTraces)
#'
#' # ProteinTraces annotation
#' head(protTraces$trace_annotation)
#' # The protein_id column from the peptide traces object becomes the new id column of the protein traces object.
#' # The last 2 columns indicate how many peptides could be observed and which were summed for quantification.
#'

proteinQuantification <- function(traces,
                                  topN = 2,
                                  keep_less = FALSE,
                                  rm_decoys = TRUE,
                                  use_sibPepCorr = FALSE,
                                  use_repPepCorr = FALSE,
                                  full_intersect_only = FALSE,
                                  quantLevel = "protein_id",
                                  verbose = TRUE, ...){
  UseMethod("proteinQuantification", traces)
}

#' @describeIn proteinQuantification Protein Quantification for single traces object
#' @export
proteinQuantification.traces <- function(traces,
                                  topN = 2,
                                  keep_less = FALSE,
                                  rm_decoys = TRUE,
                                  use_sibPepCorr = FALSE,
                                  use_repPepCorr = FALSE,
                                  full_intersect_only = FALSE,
                                  quantLevel = "protein_id",
                                  verbose = TRUE, ...){

  ## Check if it's a peptide level table
  .tracesTest(traces, "peptide")

  ## remove decoys
  if (rm_decoys) {
    n_decoys <- length(grep("^DECOY", traces$trace_annotation$protein_id))
    if ( n_decoys > 0){
      idx_decoys <- grep("^DECOY_",traces$trace_annotation$protein_id)
      traces$traces <- traces$traces[-idx_decoys]
      traces$trace_annotation<- traces$trace_annotation[-idx_decoys]
      message(n_decoys, " decoys removed")
    } else {
      message("no decoys contained/removed")
    }
  }

  if (use_sibPepCorr == TRUE) {
    if (! "SibPepCorr" %in% names(traces$trace_annotation)) {
      message("No SibPepCorr available. Please calculate SibPepCorr first.
              Function uses only intensity to select TopN peptides instead.")
      use_sibPepCorr = FALSE
    }
  }

  if(quantLevel == "protein_id"){
    message("Quantification based on 'protein_id'.")
  } else if (quantLevel == "proteoform_id") {
    if(quantLevel %in% names(traces$trace_annotation)){
      message("Quantification based on 'proteoform_id'.")
      traces$trace_annotation[,protein_id_original := protein_id]
      traces$trace_annotation[,protein_id := proteoform_id]
    } else {
      stop(paste0(quantLevel," not avilable in traces object.
      Run 'clusterPeptides' first."))
    }
  } else {
    message("The quantification level you prvided is not available.
    Quantification based on 'protein_id'.")
  }

  ## Extract wide table for calculations
  if (use_sibPepCorr == TRUE) {
    peptideTracesTable <- data.table(protein_id = traces$trace_annotation$protein_id,
                                     peptide_id = traces$trace_annotation$id,
                                     SibPepCorr = traces$trace_annotation$SibPepCorr,
                                     subset(traces$traces, select =-id))
    # Calculations in long format - sum the topN peptides per protein
    peptideTracesLong <- melt(peptideTracesTable,
                              id.vars = c("protein_id", "peptide_id", "SibPepCorr"),
                              variable.name = "fraction_number",
                              value.name = "intensity")
  } else {
    peptideTracesTable <- data.table(protein_id = traces$trace_annotation$protein_id,
                                     peptide_id = traces$trace_annotation$id,
                                     subset(traces$traces, select =-id))
    # Calculations in long format - sum the topN peptides per protein
    peptideTracesLong <- melt(peptideTracesTable,
                              id.vars = c("protein_id", "peptide_id"),
                              variable.name = "fraction_number",
                              value.name = "intensity")
  }


  peptideTracesLong[, intensity:=as.numeric(intensity)]
  peptideTracesLong[, peptide_intensity:=sum(intensity), peptide_id]
  peptideTracesLong[, n_peptides:=length(unique(peptide_id)), protein_id]
  ## the ties.method makes sure how to deal with peptides of identical intensity: "first" keeps the order of occurence
  peptideTracesLong[, peptide_intensity_rank:=rank(-peptide_intensity[1:n_peptides[1]],ties.method = "first"), protein_id]
  if (use_sibPepCorr == TRUE) {
    peptideTracesLong[, peptide_SibPepCorr_rank:=rank(-SibPepCorr[1:n_peptides[1]],ties.method = "first"), protein_id]
    peptideTracesLong[, rank_sum := peptide_intensity_rank+peptide_SibPepCorr_rank]
    peptideTracesLong[, peptide_rank:= rank(rank_sum[1:n_peptides[1]],ties.method = "first"), protein_id]
  } else {
    peptideTracesLong[, peptide_rank:= peptide_intensity_rank]
  }
  peptideTracesLong <- peptideTracesLong[peptide_rank <= topN]
  ## collect information which peptides were used for quantification
  peptideTracesLong[, quant_peptides_used:=paste(unique(peptide_id), collapse = ","), protein_id]

  ## if not wanted, kick out those with less than N peptides
  if(!keep_less){
    peptideTracesLong <- peptideTracesLong[n_peptides >= topN]
  }

  ## Sum peptides to protein level (wide) traces table
  peptideTracesTopNsumWide <- as.data.table(cast(peptideTracesLong,
                                                 protein_id ~ fraction_number,
                                                 value = "intensity",
                                                 fun.aggregate = sum))

  ## move id column to end to ensure correct quant value index
  peptideTracesTopNsumWide[, id:=protein_id]
  peptideTracesTopNsumWide <- subset(peptideTracesTopNsumWide, select=-protein_id)
  setorder(peptideTracesTopNsumWide, -id)

  ## assemble updated, protein-level trace_annotation table
  oldAnnotationPeptidelevel <- copy(traces$trace_annotation)
  #oldAnnotationPeptidelevel[,old_id:=id]
  # setnames(oldAnnotationPeptidelevel,"id","old_id")
  oldAnnotationPeptidelevel <- oldAnnotationPeptidelevel[, lapply(.SD, function(col){
    if (topN > 1) {
      #col = subset(col, select =-id)
      if(length(unique(col)) != 1){
        if(class(col) == "numeric"){
          mean(col)
        } else { #if(class(col) == "character")
          if (verbose) {
            warning(paste0("muliple entries for ",quantLevel," ", protein_id, ": ", ". Picking the first", collapse = ""))
          }
          col[1]
        }
      }else{
        col[1]
      }
    } else if (topN == 1) {
      col[old_id %in% peptideTracesTopNsumWide$id]
    }
  }),by = protein_id]
  #oldAnnotationPeptidelevel[,old_id:=NULL]

  # if ("SibPepCorr" %in% names(oldAnnotationPeptidelevel)){
  #   oldAnnotationPeptidelevel[, SibPepCorr_protein_mean:=mean(SibPepCorr), protein_id]
  #   ## removal of peptide level SibPepCorr necessary after merge.data.table update 2016-11
  #   oldAnnotationPeptidelevel = unique(subset(oldAnnotationPeptidelevel, select =-SibPepCorr))
  # }
  # if ("RepPepCorr" %in% names(oldAnnotationPeptidelevel)){
  #   oldAnnotationPeptidelevel[, RepPepCorr_protein_mean:=mean(RepPepCorr), protein_id]
  #   ## removal of peptide level RepPepCorr necessary after merge.data.table update 2016-11
  #   oldAnnotationPeptidelevel = unique(subset(oldAnnotationPeptidelevel, select =-RepPepCorr))
  # }
  # if ("meanSibPepCorr" %in% names(oldAnnotationPeptidelevel)){
  #   oldAnnotationPeptidelevel[, meanRepPepCorr_protein_mean:=mean(meanSibPepCorr), protein_id]
  #   ## removal of peptide level RepPepCorr necessary after merge.data.table update 2016-11
  #   oldAnnotationPeptidelevel = unique(subset(oldAnnotationPeptidelevel, select =-meanSibPepCorr))
  # }
  # if ("meanRepPepCorr" %in% names(oldAnnotationPeptidelevel)){
  #   oldAnnotationPeptidelevel[, meanRepPepCorr_protein_mean:=mean(meanRepPepCorr), protein_id]
  #   ## removal of peptide level RepPepCorr necessary after merge.data.table update 2016-11
  #   oldAnnotationPeptidelevel = unique(subset(oldAnnotationPeptidelevel, select =-meanRepPepCorr))
  # }
  new_annotation = unique(peptideTracesLong[,.(protein_id, n_peptides, quant_peptides_used)])
  if (quantLevel == "proteoform_id") {
    setnames(new_annotation, "n_peptides", "n_peptides_perProteoform")
  }
  common_cols <- intersect(names(oldAnnotationPeptidelevel),names(new_annotation))
  peptideTracesTopNsumWideAnnotation <- merge(oldAnnotationPeptidelevel,
                                              new_annotation,
                                              by = common_cols, all.x = FALSE, all.y = TRUE)

  peptideTracesTopNsumWideAnnotation <- unique(peptideTracesTopNsumWideAnnotation)
  peptideTracesTopNsumWideAnnotation[, id:=protein_id]
  ## Order traces alphabetically according to their id and
  setorder(peptideTracesTopNsumWideAnnotation, -id)
  setcolorder(peptideTracesTopNsumWideAnnotation, c("id", names(peptideTracesTopNsumWideAnnotation[,!"id", with = F])))
  ## assemble result protein level traces object
  traces$traces <- peptideTracesTopNsumWide
  traces$trace_annotation <- peptideTracesTopNsumWideAnnotation

  if ("genomic_coord" %in% names(traces)) {
    message("genomic_coord information only available on peptide level and will
    be removed at this point.")
    traces$genomic_coord <- NULL
  }

  ## Test and return
  if(quantLevel == "protein_id"){
    traces$trace_type <- "protein"
    .tracesTest(traces, type = "protein")
  } else if (quantLevel == "proteoform_id") {
    traces$trace_annotation[,proteoform_id := protein_id]
    traces$trace_annotation[,protein_id := protein_id_original]
    traces$trace_annotation[,protein_id_original := NULL]
    traces$trace_type <- "proteoform"
    .tracesTest(traces, type = "proteoform")
  }
  return(traces)
}


#' @describeIn proteinQuantification Protein Quantification multiple traces ojects.
#' Uses only peptides present in all traces objects for quantification.
#' @export

proteinQuantification.tracesList <- function(traces,
                                         topN = 2,
                                         keep_less = FALSE,
                                         rm_decoys = TRUE,
                                         use_sibPepCorr = FALSE,
                                         use_repPepCorr = FALSE,
                                         full_intersect_only = FALSE,
                                         quantLevel = "protein_id",
                                         verbose = TRUE, ...){
  .tracesListTest(traces, type = "peptide")

  if (full_intersect_only == TRUE) {
    intersection_peptides <- .intersect2(lapply(traces, function(x) x$traces$id))
    traces_subs <- subset(traces, trace_subset_ids = intersection_peptides)
  } else {
    traces_subs <- traces
  }

  if (quantLevel != "protein_id") {
    use_sibPepCorr = FALSE
    use_repPepCorr = FALSE
    message(paste0("Using ",quantLevel," as quantLevel doesn't support the use of
    of sibPepCorr or repPepCorr for peptide selection. Setting both options to FALSE."))
  }

  if (topN > 100) {
    traces_selected <- traces_subs
  } else {
    traces_integrated <- integrateTraceIntensities(traces_subs, aggr_corr_fun = "sum")
    if("sumSibPepCorr" %in% names(traces_integrated$trace_annotation)) {
      if("sumRepPepCorr" %in% names(traces_integrated$trace_annotation)) {
        peptideTracesTable <- data.table(protein_id = traces_integrated$trace_annotation$protein_id,
                                         peptide_id = traces_integrated$trace_annotation$id,
                                         SibPepCorr = round(traces_integrated$trace_annotation$sumSibPepCorr,digits=1),
                                         RepPepCorr = round(traces_integrated$trace_annotation$sumRepPepCorr,digits=1),
                                         subset(traces_integrated$traces, select =-id))
      } else if ("meanRepPepCorr" %in% names(traces_integrated$trace_annotation)){
        peptideTracesTable <- data.table(protein_id = traces_integrated$trace_annotation$protein_id,
                                         peptide_id = traces_integrated$trace_annotation$id,
                                         SibPepCorr = round(traces_integrated$trace_annotation$sumSibPepCorr,digits=1),
                                         RepPepCorr = round(traces_integrated$trace_annotation$meanRepPepCorr,digits=1),
                                         subset(traces_integrated$traces, select =-id))
      } else {
        traces_integrated$trace_annotation$sumRepPepCorr = 1
        peptideTracesTable <- data.table(protein_id = traces_integrated$trace_annotation$protein_id,
                                         peptide_id = traces_integrated$trace_annotation$id,
                                         SibPepCorr = round(traces_integrated$trace_annotation$sumSibPepCorr,digits=1),
                                         RepPepCorr = round(traces_integrated$trace_annotation$sumRepPepCorr,digits=1),
                                         subset(traces_integrated$traces, select =-id))
      }
    } else if ("meanSibPepCorr" %in% names(traces_integrated$trace_annotation)) {
      if("sumRepPepCorr" %in% names(traces_integrated$trace_annotation)) {
        peptideTracesTable <- data.table(protein_id = traces_integrated$trace_annotation$protein_id,
                                         peptide_id = traces_integrated$trace_annotation$id,
                                         SibPepCorr = round(traces_integrated$trace_annotation$meanSibPepCorr,digits=1),
                                         RepPepCorr = round(traces_integrated$trace_annotation$sumRepPepCorr,digits=1),
                                         subset(traces_integrated$traces, select =-id))
      } else if ("meanRepPepCorr" %in% names(traces_integrated$trace_annotation)){
        peptideTracesTable <- data.table(protein_id = traces_integrated$trace_annotation$protein_id,
                                         peptide_id = traces_integrated$trace_annotation$id,
                                         SibPepCorr = round(traces_integrated$trace_annotation$meanSibPepCorr,digits=1),
                                         RepPepCorr = round(traces_integrated$trace_annotation$meanRepPepCorr,digits=1),
                                         subset(traces_integrated$traces, select =-id))
       } else {
         traces_integrated$trace_annotation$sumRepPepCorr = 1
         peptideTracesTable <- data.table(protein_id = traces_integrated$trace_annotation$protein_id,
                                          peptide_id = traces_integrated$trace_annotation$id,
                                          SibPepCorr = round(traces_integrated$trace_annotation$meanSibPepCorr,digits=1),
                                          RepPepCorr = round(traces_integrated$trace_annotation$sumRepPepCorr,digits=1),
                                          subset(traces_integrated$traces, select =-id))
       }
    } else {
      traces_integrated$trace_annotation$sumSibPepCorr = 1
      if("sumRepPepCorr" %in% names(traces_integrated$trace_annotation)) {
        peptideTracesTable <- data.table(protein_id = traces_integrated$trace_annotation$protein_id,
                                         peptide_id = traces_integrated$trace_annotation$id,
                                         SibPepCorr = round(traces_integrated$trace_annotation$sumSibPepCorr,digits=1),
                                         RepPepCorr = round(traces_integrated$trace_annotation$sumRepPepCorr,digits=1),
                                         subset(traces_integrated$traces, select =-id))
      } else if ("meanRepPepCorr" %in% names(traces_integrated$trace_annotation)){
        peptideTracesTable <- data.table(protein_id = traces_integrated$trace_annotation$protein_id,
                                         peptide_id = traces_integrated$trace_annotation$id,
                                         SibPepCorr = round(traces_integrated$trace_annotation$sumSibPepCorr,digits=1),
                                         RepPepCorr = round(traces_integrated$trace_annotation$meanRepPepCorr,digits=1),
                                         subset(traces_integrated$traces, select =-id))
       } else {
         traces_integrated$trace_annotation$sumRepPepCorr = 1
         peptideTracesTable <- data.table(protein_id = traces_integrated$trace_annotation$protein_id,
                                          peptide_id = traces_integrated$trace_annotation$id,
                                          SibPepCorr = round(traces_integrated$trace_annotation$sumSibPepCorr,digits=1),
                                          RepPepCorr = round(traces_integrated$trace_annotation$sumRepPepCorr,digits=1),
                                          subset(traces_integrated$traces, select =-id))
       }
    }
    # Calculations in long format - sum the topN peptides per protein
    peptideTracesLong <- melt(peptideTracesTable,
                              id.vars = c("protein_id", "peptide_id", "SibPepCorr", "RepPepCorr"),
                              variable.name = "fraction_number",
                              value.name = "intensity")
    peptideTracesLong[, intensity:=as.numeric(intensity)]
    peptideTracesLong[, peptide_intensity:=sum(intensity), peptide_id]
    peptideTracesLong[, n_peptides:=length(unique(peptide_id)), protein_id]
    ## the ties.method makes sure how to deal with peptides of identical intensity: "first" keeps the order of occurence
    peptideTracesLong[, peptide_intensity_rank:=rank(-peptide_intensity[1:n_peptides[1]],ties.method = "first"), protein_id]
    if ((use_sibPepCorr == TRUE) & (use_repPepCorr == TRUE)) {
      peptideTracesLong[, peptide_SibPepCorr_rank:=.narank(-SibPepCorr[1:n_peptides[1]],ties.method = "min",na.last="keep"), protein_id]
      peptideTracesLong[, peptide_RepPepCorr_rank:=.narank(-RepPepCorr[1:n_peptides[1]],ties.method = "min",na.last="keep"), protein_id]
      peptideTracesLong[, rank_sum := peptide_intensity_rank+peptide_SibPepCorr_rank+peptide_RepPepCorr_rank]
      peptideTracesLong[, peptide_rank:= rank(rank_sum[1:n_peptides[1]],ties.method = "first"), protein_id]
    } else if ((use_sibPepCorr == TRUE) & (use_repPepCorr == FALSE)) {
      peptideTracesLong[, peptide_SibPepCorr_rank:=.narank(-SibPepCorr[1:n_peptides[1]],ties.method = "min",na.last="keep"), protein_id]
      peptideTracesLong[, rank_sum := peptide_intensity_rank+peptide_SibPepCorr_rank]
      peptideTracesLong[, peptide_rank:= rank(rank_sum[1:n_peptides[1]],ties.method = "first"), protein_id]
    } else if ((use_sibPepCorr == FALSE) & (use_repPepCorr == TRUE)) {
      peptideTracesLong[, peptide_RepPepCorr_rank:=.narank(-RepPepCorr[1:n_peptides[1]],ties.method = "min",na.last="keep"), protein_id]
      peptideTracesLong[, rank_sum := peptide_intensity_rank+peptide_RepPepCorr_rank]
      peptideTracesLong[, peptide_rank:= rank(rank_sum[1:n_peptides[1]],ties.method = "first"), protein_id]
    } else {
      peptideTracesLong[, peptide_rank:= peptide_intensity_rank]
    }
    peptideTracesLong <- peptideTracesLong[peptide_rank <= topN]
    selectedPeptides <- unique(peptideTracesLong$peptide_id)
    traces_selected <- subset(traces_subs, trace_subset_ids = selectedPeptides)
  }
  # traces_subs <- lapply(traces_subs, function(x){
  #   x$trace_annotation[,detectedIn := NULL][detected_in := NULL]
  #   x
  # })
  # class(traces_subs) <- "tracesList"
  res <- lapply(traces_selected, proteinQuantification.traces,
                topN = topN,
                keep_less = keep_less,
                rm_decoys = rm_decoys,
                use_sibPepCorr = use_sibPepCorr,
                use_repPepCorr = use_repPepCorr,
                full_intersect_only = full_intersect_only,
                quantLevel = quantLevel,
                verbose = verbose)
  class(res) <- "tracesList"
  .tracesListTest(res)
  return(res)
}

.intersect2 <- function(...) {
  args <- list(...)
  nargs <- length(args)
  if(nargs <= 1) {
    if(nargs == 1 && is.list(args[[1]])) {
      do.call(".intersect2", args[[1]])
    } else {
      stop("cannot evaluate intersection fewer than 2 arguments")
    }
  } else if(nargs == 2) {
    intersect(args[[1]], args[[2]])
  } else {
    intersect(args[[1]], .intersect2(args[-1]))
  }
}

.narank <- function(x,ties.method,na.last){
  r<-rank(x,ties.method = ties.method,na.last=na.last)
  r[is.na(x)]<-length(x)
  r
}
