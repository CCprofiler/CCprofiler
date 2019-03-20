#' Test differential association state of peptides or proteins between 2 conditions
#' @param featureVals data.table, a long-format table of intensity values within feature boundaries.
#' A featureVals table can be produced with \code{extractFeatureVals}.
#' @param compare_between Character string, The name of the column in the featureVals table
#' holding the Conditions between which differential expression should be tested.
#' @param plot logical if plot should be created, default is TRUE.
#' @param PDF logical if PDF should be created, default is FALSE.
#' @param name character string specifying the name of output if PDF=TRUE, default is "local_vs_global_stats".
#' @return A data.table containing the differential testing results for every peptide.
#' @export
testLocalVsGlobal <- function(featureVals,
                              compare_between = "Condition",
                              plot=TRUE,
                              PDF=TRUE,
                              name="local_vs_global_stats") {

featureVals_noImpute <- copy(featureVals) 
featureVals_noImpute[, intensity:=ifelse(imputedFraction==T,0,intensity)]

grpn = uniqueN(featureVals_noImpute[,.(id, feature_id, apex)])
pb <- txtProgressBar(min = 0, max = grpn, style = 3)
tests <- featureVals_noImpute[, {
  setTxtProgressBar(pb, .GRP)
  samples = unique(.SD[,get(compare_between)])
  qints = .SD[, .(s = sum(intensity)), by = .(get(compare_between), Replicate, global_intensity)]
  qints[, s_ratio:= ifelse(global_intensity==0, 0, s/global_intensity)]
  qints[, n := nrow(qints)]
  qints[, s_norm := (s_ratio * (n - 1) + 0.5)/n]
  qints[, unique_perCondition := length(unique(round(s_norm, digits = 3))), by="get"]
  n_unique_perCondition <- min(qints$unique_perCondition)
  min_nonZero <- min(length(qints[get==samples[1]][s > 0]$s),length(qints[get==samples[2]][s > 0]$s))
  if (length(unique(design_matrix$Replicate)) > 1) {
    if((n_unique_perCondition > 1) & (min_nonZero > 1)) {
      beta_model = betareg(qints$s_norm ~ qints$get)
      beta_stat = lrtest(beta_model)
      feature_mass_fraction_betaPval = beta_stat$`Pr(>Chisq)`[2]
    } else {
      feature_mass_fraction_betaPval = 2
    }
  } else {
    feature_mass_fraction_betaPval = 1
  }
  feature_mass_fraction_1 = median(qints[get==samples[1]]$s_ratio)
  feature_mass_fraction_2 = median(qints[get==samples[2]]$s_ratio)
  feature_mass_fraction_diff = feature_mass_fraction_1-feature_mass_fraction_2
  feature_mass_fraction_FC = feature_mass_fraction_1/feature_mass_fraction_2
  feature_mass_fraction_log2FC = log2(feature_mass_fraction_1/feature_mass_fraction_2)
    .(feature_mass_fraction_1 = feature_mass_fraction_1,
      feature_mass_fraction_2 = feature_mass_fraction_2,
      feature_mass_fraction_diff = feature_mass_fraction_diff,
      feature_mass_fraction_betaPval = feature_mass_fraction_betaPval,
      feature_mass_fraction_FC = feature_mass_fraction_FC,
      feature_mass_fraction_log2FC = feature_mass_fraction_log2FC
  )},
  by = .(id, feature_id, apex)]
close(pb)

tests[feature_mass_fraction_betaPval==2, feature_mass_fraction_betaPval := NA ]

if (length(unique(design_matrix$Replicate)) > 1) {
  tests$feature_mass_fraction_pBHadj <- p.adjust(tests$feature_mass_fraction_betaPval, method = "BH")
  feature_mass_fraction_pQv <- qvalue::qvalue(tests$feature_mass_fraction_betaPval, lambda = 0.4)
  tests$feature_mass_fraction_qVal <- feature_mass_fraction_pQv$qvalues
} else {
  tests$feature_mass_fraction_pBHadj <- 1
  tests$feature_mass_fraction_qVal <- 1
}

if(plot==TRUE){
  if(PDF){
    pdf(paste0(name,".pdf"), width = 3, height = 3)
  }
  hist(tests$feature_mass_fraction_1, breaks = 100)
  hist(tests$feature_mass_fraction_2, breaks = 100)
  hist(tests$feature_mass_fraction_diff, breaks = 100)
  hist(tests$feature_mass_fraction_betaPval, breaks = 100)
  hist(tests$feature_mass_fraction_pBHadj, breaks = 100)
  hist(tests$feature_mass_fraction_qVal, breaks = 100)
  if(PDF){
    dev.off()
  }
}

return(tests)
}

#' Calculate a protein-level significance from peptide level differential association tests
#' @param tests data.table, containing peptide level p-values.
#' Produced by \code{testLocalVsGlobal}.
#' A featureVals table can be produced with \code{extractFeatureVals}.
#' @return A data.table containing the differential testing results for every protein.
#' @export
aggregateLocalVsGlobal <- function(tests){
  tests[, FCpVal := (1-feature_mass_fraction_betaPval) * sign(feature_mass_fraction_diff), by=c("feature_id","apex")]
  medianPval <- tests[ ,{mPval = median(FCpVal,na.rm=T)
  .(medianPVal = 1-(mPval * sign(mPval)),
    Npeptides = .N,
    medianDiff = median(feature_mass_fraction_diff,na.rm=T)
   )},
  by = .(feature_id, apex)]
  medianPval[, pVal := pbeta(medianPVal, Npeptides/2 + 0.5, Npeptides - (Npeptides/2 + 0.5) + 1)]
  qv <- qvalue::qvalue(medianPval$pVal, lambda = 0.4)
  medianPval$qVal <- qv$qvalues
  medianPval[, pBHadj := p.adjust(pVal, method = "fdr")]
  return(medianPval)
}
