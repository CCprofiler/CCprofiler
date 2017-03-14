#' Estimate statistics based on decoy model.
#' @description Estimate statistics based on decoy model.
#' @param complex_hypothesis data.table containing complex hypotheses.
#' @param detected_complex_features data.table containing filtered complex feature results.
#' @return List with stats
#' @export
estimate_error_decoy <- function(complex_hypothesis,detected_complex_features){
  detected_complexes <- detected_complex_features$complex_id
  detected_decoys = detected_complexes[grep("d",detected_complexes)]
  detected_targets = detected_complexes[!detected_complexes %in% detected_decoys]

  all_hypotheses_tested = unique(complex_hypothesis$complex_id)
  all_decoys = all_hypotheses_tested[grep("d",all_hypotheses_tested)]
  all_targets = all_hypotheses_tested[!all_hypotheses_tested %in% all_decoys]

  FP = length(detected_decoys)
  N= length(all_decoys)
  FPR = FP/N

  TP = length(detected_targets)*(1-FPR)
  P = length(all_targets)
  TPR = TP/P

  list(FPR=FPR,TPR=TPR,TP=TP)
}

#' Estimate statistics based on manual annotation.
#' @description Estimate statistics based on manual annotation.
#' @param table data.table with complex features mapped to manual annotation.
#' @return List with stats
#' @export
estimate_errors <- function(table){
  FP_result_ids=table$results_id[which(table$manual_id_mapp == 0)]
  FN_manual_ids=table$manual_id[which(is.na(table$manual_id_mapp))]
  TP_result_ids=table$results_id[which((!is.na(table$manual_id_mapp)) & (table$manual_id_mapp != 0))]

  FP=length(FP_result_ids)
  FN=length(FN_manual_ids)
  TP=length(TP_result_ids)

  TPR=TP/(TP+FN)
  FDR=FP/(TP+FP)
  # Comparison of public peak detection algorithms for MALDI mass spectrometry data analysis
  # Chao Yang,corresponding author1 Zengyou He,1 and Weichuan Yu
  F1_score = (2*(1-FDR)*TPR)/(1-FDR+TPR)

  #high confidence features
  TP_high_result_ids = table$results_id[which((!is.na(table$manual_id_mapp)) & (table$manual_id_mapp != 0) & (table$confidence == "High"))]
  TP_high = length(TP_high_result_ids)
  FN_high_manual_ids=table$manual_id[which((is.na(table$manual_id_mapp)) & (table$confidence == "High"))]
  FN_high = length(FN_high_manual_ids)
  TPR_high = TP_high/(TP_high+FN_high)

  # compare peak apex and boundary precision
  apex_r2 = summary(lm(apex~apex.manual, data = table))$adj.r.squared
  left_pp_r2 = summary(lm(left_pp~left_pp.manual, data = table))$adj.r.squared
  right_pp_r2 = summary(lm(right_pp~right_pp.manual, data = table))$adj.r.squared
  summary_r2 = (0.5*apex_r2)+(0.25*left_pp_r2)+(0.25*right_pp_r2)

  list(
    TP=TP,
    FP=FP,
    FN=FN,
    TPR=TPR,
    FDR=FDR,
    F1_score=F1_score,
    TPR_high=TPR_high,
    summary_r2=summary_r2,
    FP_result_ids=paste(FP_result_ids,collapse=";"),
    TP_result_ids=paste(TP_result_ids,collapse=";"),
    FN_manual_ids=paste(FN_manual_ids,collapse=";")
    )
}


#' manualFeatureMapping
#' @description manualFeatureMapping
#' @param feature data.table containing filtered complex feature results.
#' @param dist_cutoff numeric maximal apex distance for features to be mapped
#' @param manual_annotation data.table containing manual annotation
#' @return data.table with mached detected and manual features
#' @export
manualFeatureMapping <- function(feature,dist_cutoff,manual_annotation){
  manual_features <- subset(manual_annotation,protein_id==feature$protein_id)
  if (nrow(manual_features) == 0) {
    match = 0
  } else {
    apex_dist <- abs(feature$apex - manual_features$apex)
    if (min(apex_dist)[1] > dist_cutoff) {
      match = 0
    } else {
      sel_min_apex_dist <- which(apex_dist==min(apex_dist))
      if (length(sel_min_apex_dist) == 1) {
        match <- manual_features$manual_id[sel_min_apex_dist]
      } else {
        manual_features <- manual_features[sel_min_apex_dist]
        left_dist <- abs(feature$left_pp - manual_features$left_pp)
        right_dist <- abs(feature$right_pp - manual_features$right_pp)
        boundary_dist <- left_dist + right_dist
        sel_min_boundary_dist <- which(boundary_dist==min(boundary_dist))
        if (length(sel_min_boundary_dist) == 1) {
          match <- manual_features$manual_id[sel_min_boundary_dist]
        } else {
          message("two manual features found")
          match <- manual_features$manual_id[sel_min_boundary_dist][1]
        }
      }
    }
  }
  match
}

#' resolveDoubleAssignments
#' @description resolveDoubleAssignments
#' @param id feature_id
#' @param mapped_features data.table with all mapped features
#' @return data.table with mached detected and manual features
#' @export
resolveDoubleAssignments <- function(id,mapped_features){
  features <- mapped_features[which(manual_id_mapp== id)]
  if (nrow(features) > 1) {
    keep <- features$results_id[which(features$peak_dist==(min(features$peak_dist)))[1]]
    unassign <- features$results_id[which(features$results_id != keep)]
    unassign
  }
}
