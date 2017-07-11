# Due to: http://stackoverflow.com/questions/24501245/data-table-throws-object-not-found-error
.datatable.aware=TRUE

#' Protein feature detection
#' @description Run the sliding window algorithm to find protein features.
#' @param traces An object of class traces (type "peptide").
#' @param corr_cutoff Numeric, the correlation value for chromatograms above which
#'        peptides are considered to be coeluting, default=0.95.
#' @param window_size Numeric, size of the window (in fractions), default=12
#' @param parallelized Logical, wether the computation should be done in parallel, default=FALSE
#' @param n_cores Integer, the number of cores to use for parallel processing
#' (only applies if parallelized is TRUE), default=1
#' @param collapse_method Method for collapsing multiple features into one feature: 
#' \itemize{
#' \item "apex_only": collapses by apex
#' \item "apex_network": collapses by apex and connected network cluster
#'}
#' Default="apex_only"
#' @param perturb_cutoff Numeric, the quantile to use in estimating the perturbation level, default="5%".
#'        Intensity values that are zero are replaced with random values that are
#'        below the specified quantile of the input values. Alternatively a
#'        cutoff value can be specified as an upper limit for perturbation values.
#'        This is nescessary for correlation calculation.
#' @param rt_height Numeric, RT cutoff for collapsing features, default is 5
#' @param smoothing_length Numeric, smoothing length of Savitzky-Golay filter, default is 7
#' @param useRandomDecoyModel Logical, wether random peptide protein associations should be used as decoy model, default = TRUE
#' @return A data.table containing protein features.
#' @export

findProteinFeatures <- function(traces,
                                corr_cutoff=0.95,
                                window_size=12,
                                parallelized=FALSE,
                                n_cores=1,
                                collapse_method="apex_only",
                                perturb_cutoff= "5%",
                                rt_height=5,
                                smoothing_length=9,
                                useRandomDecoyModel=TRUE){
  protein.peptide.association.table = traces$trace_annotation
  protein.peptide.association.table=copy(protein.peptide.association.table) # to prevent updates by reference
  setorder(protein.peptide.association.table, "protein_id")
  setnames(protein.peptide.association.table, "protein_id", "complex_id")
  setnames(protein.peptide.association.table, "id", "protein_id")
  protein.peptide.association.table[, complex_name:=complex_id]
  protein.peptide.association.table <- subset(protein.peptide.association.table,select=c("complex_id","complex_name","protein_id"))
  #protein.peptide.association.table
  if(useRandomDecoyModel){
    set.seed(123)
    random_peptides = sample(protein.peptide.association.table$protein_id) # naming is already changed for featureFinding (protein_ids == peptide_ids)
    random.association.table <- data.table(complex_id=paste0("DECOY_",protein.peptide.association.table$complex_id),
                                          complex_name=paste0("DECOY_",protein.peptide.association.table$complex_name),
                                          protein_id=paste0(random_peptides))
    protein.peptide.association.table <- rbind(protein.peptide.association.table,random.association.table)
  }

  ProteinFeatures <- findComplexFeatures(traces = traces,
                                         complex_hypothesis = protein.peptide.association.table,
                                         corr_cutoff = corr_cutoff,
                                         window_size = window_size,
                                         parallelized = parallelized,
                                         n_cores=n_cores,
                                         collapse_method=collapse_method,
                                         perturb_cutoff=perturb_cutoff,
                                         rt_height=rt_height,
                                         smoothing_length=smoothing_length)

  names(ProteinFeatures) <- gsub("complex","protein",names(ProteinFeatures))
  
  if ("protein_mw" %in% colnames(traces$trace_annotation)) {
    fun <- function(x) {strsplit(x,split=";")[[1]][1]}
    ProteinFeatures$monomer_mw <- as.numeric(unlist(lapply(ProteinFeatures$monomer_mw,fun)))
    ProteinFeatures[,monomer_mw_max := 2*monomer_mw]
    if ("molecular_weight" %in% colnames(traces$fraction_annotation)) {
      ProteinFeatures$monomer_sec <- as.numeric(unlist(lapply(ProteinFeatures$monomer_sec,fun)))
      fun <- function(X){traces$fraction_annotation$id[which.min(abs(traces$fraction_annotation$molecular_weight - X))[1]]}
      ProteinFeatures[,in_complex := FALSE]
      ProteinFeatures$in_complex[which(ProteinFeatures$apex_mw > 2*ProteinFeatures$monomer_mw)] = TRUE
      ProteinFeatures <- subset(ProteinFeatures,select=c("protein_id","protein_name","subunits_annotated","n_subunits_annotated","subunits_detected","n_subunits_detected","completeness","left_pp","right_pp","apex","apex_mw","area","peak_corr","monomer_sec","monomer_mw","in_complex"))
    } else {
      ProteinFeatures <- subset(ProteinFeatures,select=c("protein_id","protein_name","subunits_annotated","n_subunits_annotated","subunits_detected","n_subunits_detected","completeness","left_pp","right_pp","apex","area","peak_corr","monomer_mw"))
    }
  } else {
    if ("molecular_weight" %in% colnames(traces$fraction_annotation)) {
      ProteinFeatures <- subset(ProteinFeatures,select=c("protein_id","protein_name","subunits_annotated","n_subunits_annotated","subunits_detected","n_subunits_detected","completeness","left_pp","right_pp","apex","apex_mw","area","peak_corr"))
    } else {
      ProteinFeatures <- subset(ProteinFeatures,select=c("protein_id","protein_name","subunits_annotated","n_subunits_annotated","subunits_detected","n_subunits_detected","completeness","left_pp","right_pp","apex","area","peak_corr"))
    }
  }
  ProteinFeatures
}
