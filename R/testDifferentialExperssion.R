
testDifferentialExpression <- function(featureVals, 
                                       compare_between = "Condition",
                                       level = c("protein","peptide")){
  level <- match.arg(level)
  featVals <- copy(featureVals)
  setkeyv(featVals, c("feature_id", "apex", "id", "fraction"))
  message("Excluding peptides only found in one condition...")
  featureValsBoth <- filterValsByOverlap(featVals, compare_between)
  
  message("Testing peptide-level differential expression")
  
  grpn = uniqueN(featureValsBoth[,.(id, feature_id, apex)])
  pb <- txtProgressBar(min = 0, max = grpn, style = 3)
  tests <- featureValsBoth[, {
    setTxtProgressBar(pb, .GRP)
    a = t.test(formula = intensity ~ get(compare_between) , paired = T, var.equal = FALSE)
    ints = .SD[, .(s = sum(intensity)), by = get(compare_between)]$s
    int1 = ints[1]
    int2 = ints[2]
    .(pVal = a$p.value, int1 = int1, int2 = int2, meanDiff = a$estimate,
      log2FC =  log2(int1/int2),n_fractions = a$parameter + 1)},
    by = .(id, feature_id, apex)]
  close(pb)
  
  medianPval <- getFCadjustedMedian(tests)
  if(level == "peptide") return(medianPval)
  
  message("Aggregating to protein-level...")
  medianPval[, protPval := pbeta(medianPVal, (Npeptides+1)/2, 0.5*(Npeptides+1))]
  qv <- qvalue::qvalue(medianPval$protPval, lambda = 0.4)
  medianPval$QVal <- qv$qvalues
  medianPval[, pBHadj := p.adjust(protPval, method = "fdr")]
  return(medianPval)
}

getFCadjustedMedian <- function(tests){
  tests[, FCpVal := (1-pVal) * sign(log2FC)]
  medianPval <- tests[ ,{mPval = median(FCpVal)
  .(medianPVal = 1-(mPval * sign(mPval)), Npeptides = .N, medianLog2FC = median(log2FC))},
  by = .(feature_id, apex)]
  }

filterValsByOverlap <- function(featureVals, compare_between){
  # Select peptides present in both conditions
  # conditions <- unique(featureVals[,get(compare_between)])
  fv <- unique(featureVals[,.(id, feature_id, apex, get(compare_between))])
  fv$dup <- duplicated(fv[, .(id, feature_id, apex)])
  fv <- unique(fv[dup == TRUE, .(id, feature_id, apex)])
  # split <- lapply(conditions, function(cond) featureVals[get(compare_between) == cond, .(feature_id, id,apex)])
  # featureValsBoth <- featureVals[, .SD[all(sapply(conditions ,"%in%", get(compare_between)))], by = .(id, feature_id, get(compare_between))]
  featureValsBoth <- merge(featureVals, fv, by = c("id", "feature_id", "apex"))
  return(featureValsBoth)
}


normalizeVals <- function(featureVals,
                          compare_between = "Condition"){
  featureVals[, normIntensity := intensity / median(intensity[intensity > 0]), by = get(compare_between)]
  medianOfMedians <- median(featureVals[, .(median = median(intensity)), by = get(compare_between)]$median)
  featureVals[, normIntensity := normIntensity * medianOfMedians]
  return(featureVals)
}





