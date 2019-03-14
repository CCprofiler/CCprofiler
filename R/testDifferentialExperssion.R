#' Test differential expression of peptides or proteins between 2 conditions
#' @param featureVals data.table, a long-format table of intensity values within feature boundaries.
#' A featureVals table can be produced with \code{extractFeatureVals}.
#' @param compare_between Character string, The name of the column in the featureVals table
#' holding the Conditions between which differential expression should be tested.
#' @param level Character string, return tests on peptide level or aggregate to protein level.
#' Must be one of c("protein","peptide"). Defaults to "protein".
#' @param measuredOnly Logical, if only measured values should be used for differential expression testing.
#' @return A data.table containing the differential testing results for every protein/peptide.
#' @import qvalue
#' @export
testDifferentialExpression <- function(featureVals,
                                       compare_between = "Condition",
                                       level = c("protein","proteoform","peptide","complex"),
                                       measuredOnly = TRUE){
  level <- match.arg(level)
  featVals <- copy(featureVals)
  if ("complex_id" %in% names(featVals)) {
    setkeyv(featVals, c("feature_id", "complex_id", "apex", "id", "fraction"))
  } else {
    setkeyv(featVals, c("feature_id", "apex", "id", "fraction"))
  }

  message("Excluding peptides only found in one condition...")
  if (measuredOnly) {
    featVals <- subset(featVals,imputedFraction==FALSE)
    featureValsBoth <- filterValsByFractionOverlap(featVals, compare_between)
    featureValsBoth[, n_frac := .N, by=c("id", "feature_id", "apex",compare_between)]
    featureValsBoth <- subset(featureValsBoth, n_frac > 2)
  } else {
    featureValsBoth <- filterValsByOverlap(featVals, compare_between)
  }
  featureValsBoth <- getQuantTraces(featureValsBoth, compare_between)

  message("Testing peptide-level differential expression")
  if ("complex_id" %in% names(featureValsBoth)) {
    grpn = uniqueN(featureValsBoth[,.(id, feature_id, complex_id, apex)])
    pb <- txtProgressBar(min = 0, max = grpn, style = 3)
    tests <- featureValsBoth[, {
      setTxtProgressBar(pb, .GRP)
      samples = unique(.SD[,get(compare_between)])
      a = t.test(formula = intensity ~ get(compare_between) , paired = T, var.equal = FALSE)
      # qints = .SD[useForQuant == T, .(s = sum(intensity)), by = .(get(compare_between))] # this disables a lot of comparisons
      qints = .SD[, .(s = sum(intensity)), by = .(get(compare_between), Replicate)] 
      ints = .SD[imputedFraction == F, .(s = sum(intensity)), by = .(get(compare_between))] # this creates quantitative discrepancies depending on how many fractions are used
      int1 = max(0, median(ints[get==samples[1]]$s), na.rm=T)
      int2 = max(0, median(ints[get==samples[2]]$s), na.rm=T)
      qint1 = median(qints[get==samples[1]]$s)
      qint2 = median(qints[get==samples[2]]$s)
      global_ints = .SD[, .(s = unique(global_intensity)), by = .(get(compare_between), Replicate)]
      global_ints_imp = .SD[, .(s = unique(global_intensity_imputed)), by = .(get(compare_between), Replicate)]
      global_int1 = median(global_ints[get==samples[1]]$s)
      global_int2 = median(global_ints[get==samples[2]]$s)
      global_int1_imp = median(global_ints_imp[get==samples[1]]$s)
      global_int2_imp = median(global_ints_imp[get==samples[2]]$s)
      local_FC_all = log2(qints[get==samples[1]]$s/qints[get==samples[2]]$s)
      global_FC_all = log2(global_ints_imp[get==samples[1]]$s/global_ints_imp[get==samples[2]]$s)
      local_vs_global_FC_all = data.table(fc=c(local_FC_all,global_FC_all),sam=c(rep("local",length(local_FC_all)),rep("global",length(global_FC_all))))
      if (length(unique(design_matrix$Replicate)) > 2) {
        b = t.test(formula = global_ints_imp$s ~ global_ints_imp$get , paired = F, var.equal = FALSE) 
        global_pVal = b$p.value
        c = t.test(formula = local_vs_global_FC_all$fc ~ local_vs_global_FC_all$sam , paired = F, var.equal = FALSE) 
        local_vs_global_pVal = b$p.value
      } else {
        global_pVal = 1
        local_vs_global_pVal = 1
      }
      .(pVal = a$p.value, 
        int1 = int1, int2 = int2, 
        meanDiff = a$estimate,
        qint1 = qint1, qint2 = qint2, log2FC =  log2(qint1/qint2),
        n_fractions = a$parameter + 1,  Tstat = a$statistic, testOrder = paste0(samples[1],".vs.",samples[2]),
        global_int1 = global_int1, global_int2 = global_int2, global_log2FC = log2(global_int1/global_int2),
        global_int1_imp = global_int1_imp, global_int2_imp = global_int2_imp, global_log2FC_imp = log2(global_int1_imp/global_int2_imp),
        local_vs_global_log2FC = log2(qint1/qint2)-log2(global_int1/global_int2), local_vs_global_log2FC_imp = log2(qint1/qint2)-log2(global_int1_imp/global_int2_imp),
        global_pVal = global_pVal
       )},
      by = .(id, feature_id, complex_id, apex)]
    close(pb)
  } else {
    grpn = uniqueN(featureValsBoth[,.(id, feature_id, apex)])
    pb <- txtProgressBar(min = 0, max = grpn, style = 3)
    tests <- featureValsBoth[, {
      setTxtProgressBar(pb, .GRP)
      samples = unique(.SD[,get(compare_between)])
      a = t.test(formula = intensity ~ get(compare_between) , paired = T, var.equal = FALSE)
      qints = .SD[, .(s = sum(intensity)), by = .(get(compare_between), Replicate)] 
      ints = .SD[imputedFraction == F, .(s = sum(intensity)), by = .(get(compare_between))] # this creates quantitative discrepancies depending on how many fractions are used
      int1 = max(0, median(ints[get==samples[1]]$s), na.rm=T)
      int2 = max(0, median(ints[get==samples[2]]$s), na.rm=T)
      qint1 = median(qints[get==samples[1]]$s)
      qint2 = median(qints[get==samples[2]]$s)
      global_ints = .SD[, .(s = unique(global_intensity)), by = .(get(compare_between), Replicate)]
      global_ints_imp = .SD[, .(s = unique(global_intensity_imputed)), by = .(get(compare_between), Replicate)]
      global_int1 = median(global_ints[get==samples[1]]$s)
      global_int2 = median(global_ints[get==samples[2]]$s)
      global_int1_imp = median(global_ints_imp[get==samples[1]]$s)
      global_int2_imp = median(global_ints_imp[get==samples[2]]$s)
      local_FC_all = log2(qints[get==samples[1]]$s/qints[get==samples[2]]$s)
      global_FC_all = log2(global_ints_imp[get==samples[1]]$s/global_ints_imp[get==samples[2]]$s)
      local_vs_global_FC_all = data.table(fc=c(local_FC_all,global_FC_all),sam=c(rep("local",length(local_FC_all)),rep("global",length(global_FC_all))))
      if (length(unique(design_matrix$Replicate)) > 2) {
        b = t.test(formula = global_ints_imp$s ~ global_ints_imp$get , paired = F, var.equal = FALSE) 
        global_pVal = b$p.value
        c = t.test(formula = local_vs_global_FC_all$fc ~ local_vs_global_FC_all$sam , paired = F, var.equal = FALSE) 
        local_vs_global_pVal = b$p.value
      } else {
        global_pVal = 1
        local_vs_global_pVal = 1
      }
      .(pVal = a$p.value, 
        int1 = int1, int2 = int2, 
        meanDiff = a$estimate,
        qint1 = qint1, qint2 = qint2, log2FC =  log2(qint1/qint2),
        n_fractions = a$parameter + 1,  Tstat = a$statistic, testOrder = paste0(samples[1],".vs.",samples[2]),
        global_int1 = global_int1, global_int2 = global_int2, global_log2FC = log2(global_int1/global_int2),
        global_int1_imp = global_int1_imp, global_int2_imp = global_int2_imp, global_log2FC_imp = log2(global_int1_imp/global_int2_imp),
        local_vs_global_log2FC = log2(qint1/qint2)-log2(global_int1/global_int2), local_vs_global_log2FC_imp = log2(qint1/qint2)-log2(global_int1_imp/global_int2_imp),
        global_pVal = global_pVal, local_vs_global_pVal = local_vs_global_pVal
      )},
      by = .(id, feature_id, apex)]
    close(pb)
  }
  tests[is.na(log2FC) & (int1 == 0 | int2  == 0) & (meanDiff == 0)]$log2FC <- 0
  tests[is.na(log2FC) & (int1 == 0 | int2  == 0) & (meanDiff > 0)]$log2FC <- Inf
  tests[is.na(log2FC) & (int1 == 0 | int2  == 0) & (meanDiff < 0)]$log2FC <- -Inf

  if ("proteoform_id" %in% names(featVals)) {
    proteoform_ann <- unique(subset(featVals,select=c("id","proteoform_id")))
    tests <- merge(tests,proteoform_ann,by=c("id"),all.x=T,all.y=F,sort=F)
  }

  if(level == "peptide"){
    tests$pBHadj <- p.adjust(tests$pVal, method = "BH")
    pQv <- qvalue::qvalue(tests$pVal, lambda = 0.4)
    tests$qVal <- pQv$qvalues
    if (length(unique(design_matrix$Replicate)) > 2) {
      tests$global_pBHadj <- p.adjust(tests$global_pVal, method = "BH")
      global_pQv <- qvalue::qvalue(tests$global_pVal, lambda = 0.4)
      tests$global_qVal <- global_pQv$qvalues
      tests$local_vs_global_pBHadj <- p.adjust(tests$local_vs_global_pVal, method = "BH")
      local_vs_global_pQv <- qvalue::qvalue(tests$local_vs_global_pVal, lambda = 0.4)
      tests$local_vs_global_qVal <- local_vs_global_pQv$qvalues
    } else {
      tests$global_pBHadj <- 1
      tests$global_qVal <- 1
      tests$local_vs_global_pBHadj <- 1
      tests$local_vs_global_qVal <- 1
    }
    return(tests)
  } else if (level == "proteoform") {
    message("Aggregating to proteoform-level...")
    proteoformtests <- aggregatePeptideTestsToProteoform(tests)
    return(proteoformtests)
  } else if (level == "protein") {
    message("Aggregating to protein-level...")
    prottests <- aggregatePeptideTests(tests)
    return(prottests)
  } else if (level == "complex") {
    message("Aggregating to complex-level...")
    prottests <- aggregatePeptideTests(tests)
    complextests <- aggregateProteinTests(prottests)
    return(complextests)
  } else {
    stop("Specified level is not valid. Please chose between peptide, protein and complex.")
  }
}

#' Calculate a protein-level significance from peptide level differential expression tests
#' @param tests data.table, containing peptide level p-values.
#' Produced by \code{testDifferentialExpression} with option \code{level}="peptide".
#' A featureVals table can be produced with \code{extractFeatureVals}.
#' @return A data.table containing the differential testing results for every protein.
#' @export

aggregatePeptideTests <- function(tests){
  aggregated <- aggregateTests(tests,level="protein")
  return(aggregated)
}

#' Calculate a complex-level significance from protein level differential expression tests
#' @param tests data.table, containing protein level p-values.
#' Produced by \code{testDifferentialExpression} with option \code{level}="protein".
#' A featureVals table can be produced with \code{extractFeatureVals}.
#' @return A data.table containing the differential testing results for every complex.
#' @export
aggregateProteinTests <- function(tests){
  aggregated <- aggregateTests(tests,level="complex")
  return(aggregated)
}

aggregatePeptideTestsToProteoform <- function(tests){
  aggregated <- aggregateTests(tests,level="proteoform")
  return(aggregated)
}



aggregateTests <- function(tests,level){
  medianPval <- getFCadjustedMedian(tests, level=level)
  if (level == "protein") {
    medianPval[, pVal := pbeta(medianPVal, Npeptides/2 + 0.5, Npeptides - (Npeptides/2 + 0.5) + 1)]
  } else if (level == "proteoform") {
    medianPval[, pVal := pbeta(medianPVal, Npeptides/2 + 0.5, Npeptides - (Npeptides/2 + 0.5) + 1)]
  } else if (level == "complex") {
    medianPval[, pVal := pbeta(medianPVal, Nproteins/2 + 0.5, Nproteins - (Nproteins/2 + 0.5) + 1)]
  } else {
    stop("Test level must be protein, proteoform or complex.")
  }
  qv <- qvalue::qvalue(medianPval$pVal, lambda = 0.4)
  medianPval$qVal <- qv$qvalues
  medianPval[, pBHadj := p.adjust(pVal, method = "fdr")]
  return(medianPval)
}

getFCadjustedMedian <- function(tests,level){
  test <- copy(tests)
  if (level == "protein") {
    test[, FCpVal := (1-pVal) * sign(log2FC)]
    if ("complex_id" %in% names(test)) {
      medianPval <- test[ ,{mPval = median(FCpVal,na.rm=T)
      .(medianPVal = 1-(mPval * sign(mPval)),
      Npeptides = .N,
      medianLog2FC = median(log2FC,na.rm=T),
      medianTstat = median(Tstat),
      medianMeanDiff = median(meanDiff),
      qint1 = sum(qint1,na.rm=T),
      qint2 = sum(qint2,na.rm=T),
      sumLog2FC = log2(sum(qint1,na.rm=T)/sum(qint2,na.rm=T)),
      global_int1 = sum(global_int1,na.rm=T),
      global_int2 = sum(global_int2,na.rm=T),
      global_sumLog2FC = log2(sum(global_int1,na.rm=T)/sum(global_int2,na.rm=T)),
      global_int1_imp = sum(global_int1_imp,na.rm=T),
      global_int2_imp = sum(global_int2_imp,na.rm=T),
      global_sumLog2FC_imp = log2(sum(global_int1_imp,na.rm=T)/sum(global_int2_imp,na.rm=T)),
      local_vs_global_log2FC = (log2(sum(qint1,na.rm=T)/sum(qint2,na.rm=T)))-(log2(sum(global_int1,na.rm=T)/sum(global_int2,na.rm=T))),
      local_vs_global_log2FC_imp = (log2(sum(qint1,na.rm=T)/sum(qint2,na.rm=T)))-(log2(sum(global_int1_imp,na.rm=T)/sum(global_int2_imp,na.rm=T)))
      )},
      by = .(feature_id, complex_id, apex)]
    } else {
      medianPval <- test[ ,{mPval = median(FCpVal,na.rm=T)
      .(medianPVal = 1-(mPval * sign(mPval)),
      Npeptides = .N,
      medianLog2FC = median(log2FC,na.rm=T),
      medianTstat = median(Tstat),
      medianMeanDiff = median(meanDiff),
      qint1 = sum(qint1,na.rm=T),
      qint2 = sum(qint2,na.rm=T),
      sumLog2FC = log2(sum(qint1,na.rm=T)/sum(qint2,na.rm=T)),
      global_int1 = sum(global_int1,na.rm=T),
      global_int2 = sum(global_int2,na.rm=T),
      global_sumLog2FC = log2(sum(global_int1,na.rm=T)/sum(global_int2,na.rm=T)),
      global_int1_imp = sum(global_int1_imp,na.rm=T),
      global_int2_imp = sum(global_int2_imp,na.rm=T),
      global_sumLog2FC_imp = log2(sum(global_int1_imp,na.rm=T)/sum(global_int2_imp,na.rm=T)),
      local_vs_global_log2FC = (log2(sum(qint1,na.rm=T)/sum(qint2,na.rm=T)))-(log2(sum(global_int1,na.rm=T)/sum(global_int2,na.rm=T))),
      local_vs_global_log2FC_imp = (log2(sum(qint1,na.rm=T)/sum(qint2,na.rm=T)))-(log2(sum(global_int1_imp,na.rm=T)/sum(global_int2_imp,na.rm=T)))
      )},
      by = .(feature_id, apex)]
    }
  } else if (level == "complex") {
    test[, FCpVal := (1-pVal) * sign(sumLog2FC)]
    medianPval <- test[ ,{mPval = median(FCpVal,na.rm=T)
    .(medianPVal = 1-(mPval * sign(mPval)),
    Nproteins = .N,
    medianLog2FC = median(sumLog2FC,na.rm=T),
    medianTstat = median(medianTstat),
    medianMeanDiff = median(medianMeanDiff),
    qint1 = sum(qint1,na.rm=T),
    qint2 = sum(qint2,na.rm=T),
    sumLog2FC = log2(sum(qint1,na.rm=T)/sum(qint2,na.rm=T)),
    global_int1 = sum(global_int1,na.rm=T),
    global_int2 = sum(global_int2,na.rm=T),
    global_sumLog2FC = log2(sum(global_int1,na.rm=T)/sum(global_int2,na.rm=T)),
    global_int1_imp = sum(global_int1_imp,na.rm=T),
    global_int2_imp = sum(global_int2_imp,na.rm=T),
    global_sumLog2FC_imp = log2(sum(global_int1_imp,na.rm=T)/sum(global_int2_imp,na.rm=T)),
    local_vs_global_log2FC = (log2(sum(qint1,na.rm=T)/sum(qint2,na.rm=T)))-(log2(sum(global_int1,na.rm=T)/sum(global_int2,na.rm=T))),
    local_vs_global_log2FC_imp = (log2(sum(qint1,na.rm=T)/sum(qint2,na.rm=T)))-(log2(sum(global_int1_imp,na.rm=T)/sum(global_int2_imp,na.rm=T)))
    )},
    by = .(complex_id, apex)]
  } else if (level == "proteoform") {
    test[, FCpVal := (1-pVal) * sign(log2FC)]
    if ("complex_id" %in% names(test)) {
      medianPval <- test[ ,{mPval = median(FCpVal,na.rm=T)
      .(medianPVal = 1-(mPval * sign(mPval)),
      Npeptides = .N,
      medianLog2FC = median(log2FC,na.rm=T),
      medianTstat = median(Tstat),
      medianMeanDiff = median(meanDiff),
      qint1 = sum(qint1,na.rm=T),
      qint2 = sum(qint2,na.rm=T),
      sumLog2FC = log2(sum(qint1,na.rm=T)/sum(qint2,na.rm=T)),
      global_int1 = sum(global_int1,na.rm=T),
      global_int2 = sum(global_int2,na.rm=T),
      global_sumLog2FC = log2(sum(global_int1,na.rm=T)/sum(global_int2,na.rm=T)),
      global_int1_imp = sum(global_int1_imp,na.rm=T),
      global_int2_imp = sum(global_int2_imp,na.rm=T),
      global_sumLog2FC_imp = log2(sum(global_int1_imp,na.rm=T)/sum(global_int2_imp,na.rm=T)),
      local_vs_global_log2FC = (log2(sum(qint1,na.rm=T)/sum(qint2,na.rm=T)))-(log2(sum(global_int1,na.rm=T)/sum(global_int2,na.rm=T))),
      local_vs_global_log2FC_imp = (log2(sum(qint1,na.rm=T)/sum(qint2,na.rm=T)))-(log2(sum(global_int1_imp,na.rm=T)/sum(global_int2_imp,na.rm=T)))
      )},
      by = .(feature_id, proteoform_id, complex_id, apex)]
    } else {
      medianPval <- test[ ,{mPval = median(FCpVal,na.rm=T)
      .(medianPVal = 1-(mPval * sign(mPval)),
      Npeptides = .N,
      medianLog2FC = median(log2FC,na.rm=T),
      medianTstat = median(Tstat),
      medianMeanDiff = median(meanDiff),
      qint1 = sum(qint1,na.rm=T),
      qint2 = sum(qint2,na.rm=T),
      sumLog2FC = log2(sum(qint1,na.rm=T)/sum(qint2,na.rm=T)),
      global_int1 = sum(global_int1,na.rm=T),
      global_int2 = sum(global_int2,na.rm=T),
      global_sumLog2FC = log2(sum(global_int1,na.rm=T)/sum(global_int2,na.rm=T)),
      global_int1_imp = sum(global_int1_imp,na.rm=T),
      global_int2_imp = sum(global_int2_imp,na.rm=T),
      global_sumLog2FC_imp = log2(sum(global_int1_imp,na.rm=T)/sum(global_int2_imp,na.rm=T)),
      local_vs_global_log2FC = (log2(sum(qint1,na.rm=T)/sum(qint2,na.rm=T)))-(log2(sum(global_int1,na.rm=T)/sum(global_int2,na.rm=T))),
      local_vs_global_log2FC_imp = (log2(sum(qint1,na.rm=T)/sum(qint2,na.rm=T)))-(log2(sum(global_int1_imp,na.rm=T)/sum(global_int2_imp,na.rm=T)))
      )},
      by = .(feature_id, proteoform_id, apex)]
    }
    medianPval <- subset(medianPval,!is.na(proteoform_id))
  } else {
    stop("Test level must be protein or complex.")
  }
  medianPval[is.nan(sumLog2FC) & (medianLog2FC == "Inf")]$sumLog2FC <- Inf
  medianPval[is.nan(sumLog2FC) & (medianLog2FC == "-Inf")]$sumLog2FC <- -Inf
  return(medianPval)
  }

filterValsByOverlap <- function(featureVals, compare_between){
  # Select peptides present in both conditions
  # conditions <- unique(featureVals[,get(compare_between)])
  if ("complex_id" %in% names(featureVals)) {
    if("Replicate" %in% names(featureVals)){
      fv <- unique(featureVals[,.(id, feature_id, complex_id, apex, Replicate, get(compare_between))])
      fv$dup <- duplicated(fv[, .(id, feature_id, complex_id, apex, Replicate)])
      fv <- unique(fv[dup == TRUE, .(id, feature_id, complex_id, apex, Replicate)])
      featureValsBoth <- merge(featureVals, fv, by = c("id", "feature_id", "complex_id", "apex", "Replicate"))
    }else{
      fv <- unique(featureVals[,.(id, feature_id, complex_id, apex, get(compare_between))])
      fv$dup <- duplicated(fv[, .(id, feature_id, complex_id, apex)])
      fv <- unique(fv[dup == TRUE, .(id, feature_id, complex_id, apex)])
      featureValsBoth <- merge(featureVals, fv, by = c("id", "feature_id", "complex_id", "apex"))
    }
  } else {
    if("Replicate" %in% names(featureVals)){
      fv <- unique(featureVals[,.(id, feature_id, apex, Replicate, get(compare_between))])
      fv$dup <- duplicated(fv[, .(id, feature_id, apex, Replicate)])
      fv <- unique(fv[dup == TRUE, .(id, feature_id, apex, Replicate)])
      featureValsBoth <- merge(featureVals, fv, by = c("id", "feature_id", "apex", "Replicate"))
    }else{
      fv <- unique(featureVals[,.(id, feature_id, apex, get(compare_between))])
      fv$dup <- duplicated(fv[, .(id, feature_id, apex)])
      fv <- unique(fv[dup == TRUE, .(id, feature_id, apex)])
      featureValsBoth <- merge(featureVals, fv, by = c("id", "feature_id", "apex"))
    }
  }
  # split <- lapply(conditions, function(cond) featureVals[get(compare_between) == cond, .(feature_id, id,apex)])
  # featureValsBoth <- featureVals[, .SD[all(sapply(conditions ,"%in%", get(compare_between)))], by = .(id, feature_id, get(compare_between))]

  return(featureValsBoth)
}

filterValsByFractionOverlap <- function(featureVals, compare_between){
  # Select peptides present in both conditions
  # conditions <- unique(featureVals[,get(compare_between)])
  if ("complex_id" %in% names(featureVals)) {
    if("Replicate" %in% names(featureVals)){
      fv <- unique(featureVals[,.(id, feature_id, complex_id, apex, Replicate, fraction, get(compare_between))])
      fv$dup <- duplicated(fv[, .(id, feature_id, complex_id, apex, Replicate, fraction)])
      fv <- unique(fv[dup == TRUE, .(id, feature_id, complex_id, apex, Replicate, fraction)])
      featureValsBoth <- merge(featureVals, fv, by = c("id", "feature_id", "complex_id", "apex", "Replicate", "fraction"))
    }else{
      fv <- unique(featureVals[,.(id, feature_id, complex_id, apex, fraction, get(compare_between))])
      fv$dup <- duplicated(fv[, .(id, feature_id, complex_id, apex, fraction)])
      fv <- unique(fv[dup == TRUE, .(id, feature_id, complex_id, apex, fraction)])
      featureValsBoth <- merge(featureVals, fv, by = c("id", "feature_id", "complex_id", "apex", "fraction"))
    }
  } else {
    if("Replicate" %in% names(featureVals)){
      fv <- unique(featureVals[,.(id, feature_id, apex, Replicate, fraction, get(compare_between))])
      fv$dup <- duplicated(fv[, .(id, feature_id, apex, Replicate, fraction)])
      fv <- unique(fv[dup == TRUE, .(id, feature_id, apex, Replicate, fraction)])
      featureValsBoth <- merge(featureVals, fv, by = c("id", "feature_id", "apex", "Replicate", "fraction"))
    }else{
      fv <- unique(featureVals[,.(id, feature_id, apex, fraction, get(compare_between))])
      fv$dup <- duplicated(fv[, .(id, feature_id, apex, fraction)])
      fv <- unique(fv[dup == TRUE, .(id, feature_id, apex, fraction)])
      featureValsBoth <- merge(featureVals, fv, by = c("id", "feature_id", "apex", "fraction"))
    }
  }
  # split <- lapply(conditions, function(cond) featureVals[get(compare_between) == cond, .(feature_id, id,apex)])
  # featureValsBoth <- featureVals[, .SD[all(sapply(conditions ,"%in%", get(compare_between)))], by = .(id, feature_id, get(compare_between))]

  return(featureValsBoth)
}

getQuantTraces <- function(featureVals, compare_between){
  if ("complex_id" %in% names(featureVals)) {
    if("Replicate" %in% names(featureVals)){
      featureVals[, useForQuant := (!any(imputedFraction) & .N == 2),
                  by=.(id, feature_id, complex_id, apex, Replicate, fraction)]
    }else{
      featureVals[, useForQuant := (!any(imputedFraction) & .N == 2),
                  by=.(id, feature_id, complex_id, apex, fraction)]
    }
  } else {
    if("Replicate" %in% names(featureVals)){
      featureVals[, useForQuant := (!any(imputedFraction) & .N == 2),
                  by=.(id, feature_id, apex, Replicate, fraction)]
    }else{
      featureVals[, useForQuant := (!any(imputedFraction) & .N == 2),
                  by=.(id, feature_id, apex, fraction)]
    }
  }
  return(featureVals)
}

normalizeVals <- function(featureVals,
                          compare_between = "Condition"){
  featureVals[, normIntensity := intensity / median(intensity[intensity > 0]), by = get(compare_between)]
  medianOfMedians <- median(featureVals[, .(median = median(intensity)), by = get(compare_between)]$median)
  featureVals[, normIntensity := normIntensity * medianOfMedians]
  return(featureVals)
}

#' Volcano plot for differential expression of peptides, proteins or complexes between 2 conditions
#' @param testResults data.table, a data.table with test statistics.
#' A testResults can be produced with \code{testDifferentialExpression}.
#' @param highlight character string of feature_id that should be highlighted. Default = NULL.
#' @param FC_cutoff Numeric fold change cutoff, default is 2.
#' @param pBHadj_cutoff Numeric p-value cutoff, default is 0.01.
#' @param name character string specifying the name of output if PDF=TRUE, default is "volcanoPlot".
#' @param PDF logical if PDF should be created, default is FALSE.
#' @return plot
#' @export
plotVolcano <- function(testResults, highlight=NULL, FC_cutoff=2, pBHadj_cutoff=0.01,name="volcanoPlot", PDF=FALSE) {
  if (PDF) {
    pdf(paste0(name,".pdf"), height=4, width=4)
  }
  if ("sumLog2FC" %in% names(testResults)) {
    p <- ggplot(testResults, aes(x=sumLog2FC,y=-log10(pBHadj)))
  } else {
    p <- ggplot(testResults, aes(x=log2FC,y=-log10(pBHadj)))
  }
  p <- p +
    geom_point(size=1) +
    theme_classic() +
    geom_hline(yintercept=-log10(pBHadj_cutoff), colour="red", linetype="dashed") +
    geom_vline(xintercept=-log2(FC_cutoff), colour="red", linetype="dashed") +
    geom_vline(xintercept=log2(FC_cutoff), colour="red", linetype="dashed")
  if (! is.null(highlight)){
    if ("feature_id" %in% names(testResults)) {
      sub <- subset(testResults,feature_id %in% highlight)
      col <- "feature_id"
    } else if ("complex_id" %in% names(testResults)) {
      sub <- subset(testResults,complex_id %in% highlight)
      col <- "complex_id"
    } else if (highlight %in% testResults$protein_id) {
      sub <- subset(testResults,protein_id %in% highlight)
      col <- "protein_id"
    } else if (highlight %in% testResults$proteoform_id) {
      sub <- subset(testResults,proteoform_id %in% highlight)
      col <- "proteoform_id"
    } else {
      stop("The testResults do not have the proper format. Input should be the result from testDifferentialExpression.")
    }
    if ("sumLog2FC" %in% names(testResults)) {
      p <- p + geom_point(data=sub, aes(x=sumLog2FC,y=-log10(pBHadj)), colour="red", fill="red", size=3, shape=23) +
        geom_text_repel(data=sub, aes(label=get(col)), size=4, vjust=0, hjust=-0.1, colour="red")
    } else {
      p <- p + geom_point(data=sub, aes(x=log2FC,y=-log10(pBHadj)), colour="red", fill="red", size=3, shape=23)+
        geom_text_repel(data=sub, aes(label=get(col)), size=4, vjust=0, hjust=-0.1, colour="red")
    }
  }
  print(p)
  if (PDF) {
    dev.off()
  }
}
