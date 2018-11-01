#' Perform global differential expression analysis
#' @param traces object of class traces
#' @param design_matrix design matrix
#' @return list with peptide and protein level differential expression results.
#' If proteoform_ids are available also proteoform differential expression in
#' estimated.
#' @export

testGlobalDifferentialExpression <- function(traces,design_matrix=NULL){
  total_fractions.dt <- lapply(traces, "[[", 4)
  total_fractions <- max(unique(unlist(lapply(total_fractions.dt, function(x){x$id}))))

  ann <- lapply(traces, "[[", 3)
  ann_all <- unique(do.call(rbind,lapply(ann, function(x){subset(x,select=c("id","protein_id"))})))
  ann_all[,protein_name := protein_id]
  ann_all[,subunits_annotated := paste(id, collapse=';'), by="protein_id"]
  mock_features <- unique(ann_all[,id:=NULL])
  mock_features[,n_subunits_annotated:=.N, by="protein_id"]
  mock_features[,subunits_detected:=subunits_annotated]
  mock_features[,n_subunits_detected:=n_subunits_annotated]
  mock_features[,completeness:=1]
  mock_features[,left_pp:=1]
  mock_features[,right_pp:=total_fractions]
  mock_features[,apex:=1]
  mock_features[,apex_mw:=1]
  mock_features[,area:=1]
  mock_features[,peak_corr:=1]
  mock_features[,monomer_sec:=1]
  mock_features[,monomer_mw:=1]
  mock_features[,in_complex:=FALSE]
  mock_features[,coelution_score:=1]
  mock_features[,decoy:=1]
  mock_features[,pvalue:=1]
  mock_features[,qvalue:=1]

  mock_features_featureVals <- extractFeatureVals(traces = traces,
                                           features = mock_features,
                                           design_matrix = design_matrix,
                                           extract = "subunits_detected",
                                           imputeZero = T,
                                           verbose = F,
                                           perturb_cutoff = "5%")

  mock_features_featureValsFilled <- fillFeatureVals(featureVals = mock_features_featureVals,
                                              design_matrix = design_matrix)

  mock_features_DiffExprPep <- testDifferentialExpression(featureVals = mock_features_featureValsFilled,
                                                   compare_between = "Condition",
                                                   level = "peptide",
                                                   measuredOnly = FALSE)

  mock_features_DiffExprProtein <- aggregatePeptideTests(mock_features_DiffExprPep)

  diffProteins <- subset(mock_features_DiffExprProtein,select=c("feature_id","pVal","global_int1_imp","global_int2_imp","global_sumLog2FC_imp","pBHadj","qVal"))
  setnames(diffProteins,"feature_id","protein_id")

  if("proteoform_id" %in% names(traces[[1]]$trace_annotation)) {
    diffPeptides <- subset(mock_features_DiffExprPep,select=c("id","feature_id","proteoform_id","pVal","global_int1_imp","global_int2_imp","global_log2FC_imp","pBHadj","qVal"))
    setnames(diffPeptides,"feature_id","protein_id")
    mock_features_DiffExprProteoform <- aggregatePeptideTestsToProteoform(mock_features_DiffExprPep)
    diffProteoforms <- subset(mock_features_DiffExprProteoform,select=c("proteoform_id","feature_id","pVal","global_int1_imp","global_int2_imp","global_sumLog2FC_imp","pBHadj","qVal"))
    setnames(diffProteoforms,"feature_id","protein_id")
    return(list(diffPeptides=diffPeptides,diffProteins=diffProteins,diffProteoforms=diffProteoforms))
  } else {
    diffPeptides <- subset(mock_features_DiffExprPep,select=c("id","feature_id","pVal","global_int1_imp","global_int2_imp","global_log2FC_imp","pBHadj","qVal"))
    setnames(diffProteins,"feature_id","protein_id")
    return(list(diffPeptides=diffPeptides,diffProteins=diffProteins))
  }
}
