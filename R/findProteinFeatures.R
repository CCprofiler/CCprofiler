#' Protein feature detection
#' @description Run the sliding window algorithm to find protein features.
#' @param pepTraces An object of type \code{traces.obj}.
#' @param corr.cutoff The correlation value for chromatograms above which
#'        peptides are considered to be coeluting.
#' @param window.size Size of the window. Numeric.
#' @param parallelized If the computation should be done in parallel.
#'        Optional.
#' @param n.cores The number of cores to use for parallel processing. Optional.
#' @param perturb.cutoff The quantile to use in estimating the perturbation level.
#'        Intensity values that are zero are replaced with random values that are
#'        below the specified quantile of the input values. Alternatively a
#'        cutoff value can be specified as an upper limit for perturbation values.
#'        This is nescessary for correlation calculation.
#' @return A data.table containing several columns TODO.
#' @export

findProteinFeatures <- function(pepTraces,
                                MWSECcalibrationFunctions,
                                corr.cutoff,
                                window.size,
                                parallelized,
                                n.cores,
                                collapse_method,
                                perturb.cutoff,
                                rt_height=5,
                                smoothing_length=11){
  protein.peptide.association.table <- as.data.table(as.data.frame(pepTraces[["trace_annotation"]]))
  setorder(protein.peptide.association.table, "protein_id")
  setnames(protein.peptide.association.table, "protein_id", "complex_id")
  setnames(protein.peptide.association.table, "id", "protein_id")
  protein.peptide.association.table[, complex_name:=complex_id]
  protein.peptide.association.table <- subset(protein.peptide.association.table,select=c("complex_id","complex_name","protein_id"))
  #protein.peptide.association.table

  ProteinFeatures <- findComplexFeatures(traces.obj = pepTraces,
                                         complex.protein.assoc = protein.peptide.association.table,
                                         MWSECcalibrationFunctions = MWSECcalibrationFunctions,
                                         corr.cutoff = corr.cutoff,
                                         window.size = window.size,
                                         parallelized = parallelized,
                                         n.cores=n.cores,
                                         collapse_method=collapse_method,
                                         perturb.cutoff=perturb.cutoff,
                                         rt_height=rt_height,
                                         smoothing_length=smoothing_length)
  protRes <- resultsToTable(ProteinFeatures)

  fun <- function(x) {strsplit(x,split=";")[[1]][1]}
  protRes$monomer_mw <- as.numeric(unlist(lapply(protRes$monomer_mw,fun)))
  protRes$monomer_sec <- as.numeric(unlist(lapply(protRes$monomer_sec,fun)))
  names(protRes) <- gsub("complex","protein",names(protRes))

  protRes <- subset(protRes,select=c("protein_id","protein_name","subunits_annotated","n_subunits_annotated","subunits_detected","n_subunits_detected","completeness","left_pp","right_pp","apex","apex_mw","area","peak_corr","monomer_sec","monomer_mw"))
  protRes[,monomer_mw_max := 2*protRes$monomer_mw]
  protRes[,monomer_sec_max := MWSECcalibrationFunctions$MWtoSECfraction(monomer_mw_max)]
  protRes[,in_complex := FALSE]
  protRes$in_complex[protRes$apex_mw > 2*protRes$monomer_mw] = TRUE

  protRes
}
