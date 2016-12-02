#' Run SECprofiler from SECexplorer y providing ids and their type.
#' @param ids A character vector containing names or identifiers.
#' @param type A character string specifying the type of ids.
#' @return A list containing the following elements:
#'    \itemize{
#'      \item \code{conversionRes} A list of the following:
#'        \itemize{
#'           \item \code{not_defined} A character vector containing all input ids that are not defined in the current version of the package.
#'           \item \code{no_ms_signal} A chaacter vector of all input ids that do not have an MS signal in the current SEC experiment.
#'           \item \code{protein_ids}  A chaacter vector of all uniprotkb accession ids used for further analysis.
#'           \item \code{mapping_table} A data.table containing the following information:
#'           \itemize{
#'            \item \code{UNIPROTKB} Uniprotkb accession ids.
#'            \item \code{name} The combined name for the uniprotkb accession ids. This is \code{UNIPROTKB} | uniprot name.
#'            \item \code{input type} If the input type was not \code{UNIPROTKB} this contains the ids of the input type.
#'           }
#'        }
#'      \item \code{featureRes} A list of the following:
#'      \itemize{
#'        \item \code{traces} An object of type \code{traces.obj}. Filtered for input ids.
#'        \item \code{features} An object of type \code{complexFeaturesAnnotated} that is a list containing the following:
#'        \itemize{
#'          \item \code{feature} data.table containing complex feature candidates in the following format:
#'           \itemize{
#'           \item \code{complex_id} The complex_id of the query complex.
#'           \item \code{complex_name} The complex_name of the query complex.
#'           \item \code{subunits_annotated} The subunits (protein_ids) annotated for the complex separated by semi-colons.
#'           \item \code{n_subunits_annotated} The number of subunits (protein_ids) annotated for the complex separated by semi-colons.
#'           \item \code{subunits_with_signal} The subunits (protein_ids) with an MS/MS signal for the complex separated by semi-colons.
#'           \item \code{n_subunits_with_signal} The number of subunits (protein_ids) with an MS/MS signal for the complex separated by semi-colons.
#'           \item \code{subunits_detected} The subunits (protein_ids) detceted in the feature for the complex separated by semi-colons.
#'           \item \code{n_subunits_detected} The number of subunits (protein_ids) detceted in the feature for the complex separated by semi-colons.
#'           \item \code{completeness} The complex completeness as defined by \code{n_subunits_detected} divided by \code{n_subunits_annotated}.
#'           \item \code{left_sw} The left boundary of the sliding-window feature.
#'           \item \code{right_sw} The right boundary of the sliding-window feature.
#'           \item \code{sw_score} The intra-sliding-window-feature correlation.
#'           \item \code{left_pp} The left boundary of the selected peak by the peak-picker.
#'           \item \code{right_pp} The right boundary of the selected peak by the peak-picker.
#'           \item \code{apax} The apex of the selected peak by the peak-picker.
#'           \item \code{area} The area (entire complex) of the selected peak by the peak-picker.
#'           \item \code{total_intensity} The intensity of all protein_ids of the feature separated by semi-colons.
#'           \item \code{intensity_ratio} The intensity ratio of all protein_ids of the feature separated by semi-colons.
#'           \item \code{stoichiometry_estimated} The rounded \code{intensity_ratio} of all protein_ids of the feature separated by semi-colons.
#'           \item \code{monomer_mw} The monomer molecular weights of all protein_ids of the feature separated by semi-colons.
#'           \item \code{monomer_sec} The monomer sec fraction of all protein_ids of the feature separated by semi-colons.
#'           \item \code{complex_mw_estimated} The complex molecular weight as expected fro the \code{stoichiometry_estimated}.
#'           \item \code{complex_sec_estimated} The complex sec fraction as expected fro the \code{stoichiometry_estimated}.
#'           \item \code{sec_diff} Difference between \code{complex_sec_estimated} and \code{apax} of the feature.
#'          }
#'        }
#'        }
#'        }
#' @export
runSECexplorer <- function(ids,type){
  traces.obj <- traces.obj
  conversionRes <- convertIDs(ids,type,traces.obj)
  featureRes <- SECexplorer_processing(conversionRes$protein_ids,traces.obj)
  list(conversionRes = conversionRes,featureRes = featureRes)
}

#protein_ids = c("P20248","P24941","P18858","P09874","P09884","P28340","Q07864","P35251","P35250","P27694","P15927","P35244","P11387")
#protein_ids = c("P61201","Q9UNS2","Q9BT78","Q92905","Q7L5N1","Q99627","Q13098","Q9UBW8","P61201","Q9UNS2","Q9BT78","Q92905","Q7L5N1","Q99627","Q13098","Q9H9Q2")

#' SECexplorer processing.
#' @param protein_ids A character vector containing \code{UNIPROTKB} accession ids.
#' @param traces.obj An object of type \code{traces.obj}.
#' @return A list containing the following:
#'      \itemize{
#'        \item \code{traces} An object of type \code{traces.obj}. Filtered for input ids.
#'        \item \code{features} An object of type \code{complexFeaturesAnnotated} that is a list containing the following:
#'        \itemize{
#'          \item \code{feature} data.table containing complex feature candidates in the following format:
#'           \itemize{
#'           \item \code{complex_id} The complex_id of the query complex.
#'           \item \code{complex_name} The complex_name of the query complex.
#'           \item \code{subunits_annotated} The subunits (protein_ids) annotated for the complex separated by semi-colons.
#'           \item \code{n_subunits_annotated} The number of subunits (protein_ids) annotated for the complex separated by semi-colons.
#'           \item \code{subunits_with_signal} The subunits (protein_ids) with an MS/MS signal for the complex separated by semi-colons.
#'           \item \code{n_subunits_with_signal} The number of subunits (protein_ids) with an MS/MS signal for the complex separated by semi-colons.
#'           \item \code{subunits_detected} The subunits (protein_ids) detceted in the feature for the complex separated by semi-colons.
#'           \item \code{n_subunits_detected} The number of subunits (protein_ids) detceted in the feature for the complex separated by semi-colons.
#'           \item \code{completeness} The complex completeness as defined by \code{n_subunits_detected} divided by \code{n_subunits_annotated}.
#'           \item \code{left_sw} The left boundary of the sliding-window feature.
#'           \item \code{right_sw} The right boundary of the sliding-window feature.
#'           \item \code{sw_score} The intra-sliding-window-feature correlation.
#'           \item \code{left_pp} The left boundary of the selected peak by the peak-picker.
#'           \item \code{right_pp} The right boundary of the selected peak by the peak-picker.
#'           \item \code{apax} The apex of the selected peak by the peak-picker.
#'           \item \code{area} The area (entire complex) of the selected peak by the peak-picker.
#'           \item \code{total_intensity} The intensity of all protein_ids of the feature separated by semi-colons.
#'           \item \code{intensity_ratio} The intensity ratio of all protein_ids of the feature separated by semi-colons.
#'           \item \code{stoichiometry_estimated} The rounded \code{intensity_ratio} of all protein_ids of the feature separated by semi-colons.
#'           \item \code{monomer_mw} The monomer molecular weights of all protein_ids of the feature separated by semi-colons.
#'           \item \code{monomer_sec} The monomer sec fraction of all protein_ids of the feature separated by semi-colons.
#'           \item \code{complex_mw_estimated} The complex molecular weight as expected fro the \code{stoichiometry_estimated}.
#'           \item \code{complex_sec_estimated} The complex sec fraction as expected fro the \code{stoichiometry_estimated}.
#'           \item \code{sec_diff} Difference between \code{complex_sec_estimated} and \code{apax} of the feature.
#'          }
#'        }
#'        }
SECexplorer_processing <- function(protein_ids,traces.obj){
  #@TODO ID mapping?
  traces.obj <- subset(traces.obj,trace_ids = protein_ids)
  complex.table <- data.table(complex_id=rep(1,length(protein_ids)),complex_name=rep("1",length(protein_ids)),protein_id=protein_ids)
  std_weights_kDa = c(1398, 699, 300, 150, 44, 17)
  std_elu_fractions = c(19, 29, 37, 46, 54.5, 61)
  calibration = calibrateSECMW(std_weights_kDa = std_weights_kDa, std_elu_fractions = std_elu_fractions,plot=TRUE,PDF=FALSE)
  swf <- findComplexFeatures(traces.obj = traces.obj,
                             complex.protein.assoc = complex.table,
                             MWSECcalibrationFunctions=calibration,
                             corr.cutoff=0.95,
                             window.size=15,
                             parallelized=FALSE,
                             perturb.cutoff = "5%",
                             collapse_method="apex_only",
                             rt_height=5,
                             smoothing_length=11)
  res = resultsToTable(swf)
  model = lm(log(std_weights_kDa) ~ std_elu_fractions)
  model_coeffitients = model$coefficients
  names(model_coeffitients) = c("intercept","slope")
  list(traces=traces.obj,features=res,calibration=model_coeffitients)
}

#TEST
#type = "UNIPROTKB"
#ids = c("P20248","BLA","Q9P1F3","Q12979","P09884","P28340","Q07864")
#type = "UNIPROTKB_name"
#ids=c("CCNA2_HUMAN","BLA","PAL4E_HUMAN","S23IP_HUMAN")
#type="Gene_names"
#ids=c("ABR","ABRACL","BLA","C6orf115")
#traces.obj <- traces.obj
#convertIDs(ids,type,traces.obj)

#' Convert any input ids to uniprotkb accession ids.
#' @param ids A character vector containing names or identifiers.
#' @param type A character string specifying the type of ids.
#' @param traces.obj An object of type \code{traces.obj}.
#' @return A data.table containing the following information:
#'        \itemize{
#'           \item \code{not_defined} A character vector containing all input ids that are not defined in the current version of the package.
#'           \item \code{no_ms_signal} A chaacter vector of all input ids that do not have an MS signal in the current SEC experiment.
#'           \item \code{protein_ids}  A chaacter vector of all uniprotkb accession ids used for further analysis.
#'           \item \code{mapping_table} A data.table containing the following information:
#'           \itemize{
#'            \item \code{UNIPROTKB} Uniprotkb accession ids.
#'            \item \code{name} The combined name for the uniprotkb accession ids. This is \code{UNIPROTKB} | uniprot name.
#'            \item \code{input type} If the input type was not \code{UNIPROTKB} this contains the ids of the input type.
#'           }
#'        }
convertIDs <- function(ids=ids,type="UNIPROTKB",traces.obj=traces.obj){
  human_annotated_more <- uniprot_reviewed_9606_10082016_annotated
  if(type=="UNIPROTKB"){
    # which input names are not defined in our mapping table as UNIPROTKB
    not_defined <- ids[which(! ids %in% human_annotated_more$UNIPROTKB)]
    # remaining protein uniprot ids
    protein_ids=unique(human_annotated_more$UNIPROTKB[which(human_annotated_more$UNIPROTKB %in% ids)])
    # which input ids do not have an MS signal in the traces object
    no_ms_signal <- protein_ids[which(! protein_ids %in% traces.obj$traces$id)]
    # final list of uniprot accessions used for the analysis
    protein_ids <- protein_ids[which(protein_ids %in% traces.obj$traces$id)]
  } else if (type=="UNIPROTKB_name") {
    # which input names are not defined in our mapping table as UNIPROTKB_name
    not_defined <- ids[which(! ids %in% human_annotated_more$UNIPROTKB_name)]
    # remaining protein uniprot ids
    protein_ids=unique(human_annotated_more$UNIPROTKB[which(human_annotated_more$UNIPROTKB_name %in% ids)])
    # which input ids do not have an MS signal in the traces object
    no_ms_signal <- protein_ids[which(! protein_ids %in% traces.obj$traces$id)]
    # convert back to the input UNIPROTKB_name
    no_ms_signal <- unique(human_annotated_more$UNIPROTKB_name[which((human_annotated_more$UNIPROTKB %in% no_ms_signal) & (human_annotated_more$UNIPROTKB_name %in% ids))])
    # final list of uniprot accessions used for the analysis
    protein_ids <- protein_ids[which(protein_ids %in% traces.obj$traces$id)]
  } else if (type=="Gene_names") {
    # which input names are not defined in our mapping table as Gene_names
    not_defined <- ids[which(! ids %in% human_annotated_more$Gene_names)]
    # remaining protein uniprot ids
    protein_ids=unique(human_annotated_more$UNIPROTKB[which(human_annotated_more$Gene_names %in% ids)])
    # which input ids do not have an MS signal in the traces object
    no_ms_signal <- protein_ids[which(! protein_ids %in% traces.obj$traces$id)]
    # convert back to the input Gene_names
    no_ms_signal <- unique(human_annotated_more$Gene_names[which((human_annotated_more$UNIPROTKB %in% no_ms_signal) & (human_annotated_more$Gene_names %in% ids))])
    # final list of uniprot accessions used for the analysis
    protein_ids <- protein_ids[which(protein_ids %in% traces.obj$traces$id)]
  } else if (type=="ENSEMBL") {
    # which input names are not defined in our mapping table as ENSEMBL
    not_defined <- ids[which(! ids %in% human_annotated_more$ENSEMBL)]
    # remaining protein uniprot ids
    protein_ids=unique(human_annotated_more$UNIPROTKB[which(human_annotated_more$ENSEMBL %in% ids)])
    # which input ids do not have an MS signal in the traces object
    no_ms_signal <- protein_ids[which(! protein_ids %in% traces.obj$traces$id)]
    # convert back to the input ENSEMBL
    no_ms_signal <- unique(human_annotated_more$ENSEMBL[which((human_annotated_more$UNIPROTKB %in% no_ms_signal) & (human_annotated_more$ENSEMBL %in% ids))])
    # final list of uniprot accessions used for the analysis
    protein_ids <- protein_ids[which(protein_ids %in% traces.obj$traces$id)]
  } else if (type=="HGNC") {
    # which input names are not defined in our mapping table as HGNC
    not_defined <- ids[which(! ids %in% human_annotated_more$HGNC)]
    # remaining protein uniprot ids
    protein_ids=unique(human_annotated_more$UNIPROTKB[which(human_annotated_more$HGNC %in% ids)])
    # which input ids do not have an MS signal in the traces object
    no_ms_signal <- protein_ids[which(! protein_ids %in% traces.obj$traces$id)]
    # convert back to the input HGNC
    no_ms_signal <- unique(human_annotated_more$HGNC[which((human_annotated_more$UNIPROTKB %in% no_ms_signal) & (human_annotated_more$HGNC %in% ids))])
    # final list of uniprot accessions used for the analysis
    protein_ids <- protein_ids[which(protein_ids %in% traces.obj$traces$id)]
  }
  if(type=="UNIPROTKB"){
    mapping_table <- unique(subset(human_annotated_more,select=c("UNIPROTKB","name"),UNIPROTKB %in% protein_ids))
  } else {
    mapping_table <- unique(subset(human_annotated_more,select=c("UNIPROTKB","name",type),UNIPROTKB %in% protein_ids))
  }
  return(list(not_defined=not_defined, no_ms_signal=no_ms_signal, protein_ids=protein_ids,mapping_table=mapping_table))
}
