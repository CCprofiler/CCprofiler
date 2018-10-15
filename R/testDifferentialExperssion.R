#' Test differential expression of peptides or proteins between 2 conditions
#' @param featureVals data.table, a long-format table of intensity values within feature boundaries.
#' A featureVals table can be produced with \code{extractFeatureVals}.
#' @param compare_between Character string, The name of the column in the featureVals table
#' holding the Conditions between which differential expression should be tested.
#' @param level Character string, return tests on peptide level or aggregate to protein level.
#' Must be one of c("protein","peptide"). Defaults to "protein".
#' @param transformation arithmetic function, the mathematical operation to transform the
#' 'intensity' values with. E.g. log, sqrt, etc. If set to NULL no transformation
#'  is performed (default).
#' @return A data.table containing the differential testing results for every protein/peptide.
#' @import qvalue
#' @export
testDifferentialExpression <- function(featureVals, 
                                       compare_between = "Condition",
                                       level = c("protein","peptide"),
                                       transformation = NULL){
  level <- match.arg(level)
  featVals <- copy(featureVals)

  setkeyv(featVals, c("feature_id", "apex", "id", "fraction"))
  message("Excluding peptides only found in one condition...")
  featureValsBoth <- filterValsByOverlap(featVals, compare_between)

  message("Summarizing peptides on feature level")
  featureValsBothSum <- unique(featureValsBoth[,.(intensity=sum(intensity),
                                                  get(compare_between)),
                                        by=.(Sample, id, feature_id, apex)])
  setnames(featureValsBothSum, "V2", compare_between)

  ## Transformation with specified function
  if (!is.null(transformation)) { #fix 1 10.5.18 MH
    if (class(transformation) == "character") {
      message(paste("Transforming intensities with", transformation, 
                    "..."))
      transformation <- get(transformation)
    }
    else {
      message(paste("Transforming intensities..."))
    }
    featureValsBothSum[, `:=`(intensityTr, transformation(intensity + 1))]
  }
  else {
    featureValsBothSum[, intensityTr:=intensity+1] # fix2 10.5.18 MH
  }
  setkeyv(featVals, c("feature_id", "apex", "id", "fraction"))
  message("Excluding peptides only found in one condition...")
  featureValsBoth <- filterValsByOverlap(featVals, compare_between)
  
  message("Testing peptide-level differential expression")
  
  grpn = uniqueN(featureValsBoth[,.(id, feature_id, apex)])
  pb <- txtProgressBar(min = 0, max = grpn, style = 3)
  tests <- featureValsBothSum[, {
    setTxtProgressBar(pb, .GRP)
    a = t.test(formula = intensityTr ~ get(compare_between), 
               paired = F, var.equal = FALSE)
    ints = .SD[, .(s = sum(intensity)), by = .(get(compare_between))]$s
    int1 = ints[1]
    int2 = ints[2]
    .(pVal = a$p.value, int1 = int1, int2 = int2, meanDiff = a$estimate, 
      FC = int1/int2, log2FC = log2(int1/int2), # Fix3: make log2FC column for downstream plot
      n_fractions = a$parameter + 1, Tstat = a$statistic)
  }, by = .(id, feature_id, apex)]
  close(pb)
  if (level == "peptide") {
    tests$pBHadj <- p.adjust(tests$pVal, method = "BH")
    pQv <- qvalue::qvalue(tests$pVal, lambda = 0.4)
    tests$QVal <- pQv$qvalues
    return(tests)
  }
  message("Aggregating to protein-level...")
  prottests <- aggregatePeptideTests(tests)
  return(prottests)
}

#' Calculate a protein-level significance from peptide level differential expression tests
#' @param tests data.table, containing peptide level p-values.
#' Produced by \code{testDifferentialExpression} with option \code{level}="peptide".
#' A featureVals table can be produced with \code{extractFeatureVals}.
#' @return A data.table containing the differential testing results for every protein.
#' @export

aggregatePeptideTests <- function(tests){
  medianPval <- getFCadjustedMedian(tests)
  medianPval[, protPval := pbeta(medianPVal, Npeptides/2 + 0.5, Npeptides - (Npeptides/2 + 0.5) + 1)]
  qv <- qvalue::qvalue(medianPval$protPval, lambda = 0.4)
  medianPval$QVal <- qv$qvalues
  medianPval[, pBHadj := p.adjust(protPval, method = "fdr")]
  return(medianPval)
  
}

getFCadjustedMedian <- function(tests){
  test <- copy(tests)
  test[, FCpVal := (1-pVal) * sign(log2(FC))]
  medianPval <- test[ ,{mPval = median(FCpVal)
  .(medianPVal = 1-(mPval * sign(mPval)), Npeptides = .N, medianFC = median(FC), medianTstat = median(Tstat))},
  by = .(feature_id, apex)]
  }

filterValsByOverlap <- function(featureVals, compare_between){
  # Select peptides present in both conditions
  # conditions <- unique(featureVals[,get(compare_between)])
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
  # split <- lapply(conditions, function(cond) featureVals[get(compare_between) == cond, .(feature_id, id,apex)])
  # featureValsBoth <- featureVals[, .SD[all(sapply(conditions ,"%in%", get(compare_between)))], by = .(id, feature_id, get(compare_between))]
  
  return(featureValsBoth)
}


normalizeVals <- function(featureVals,
                          compare_between = "Condition"){
  featureVals[, normIntensity := intensity / median(intensity[intensity > 0]), by = get(compare_between)]
  medianOfMedians <- median(featureVals[, .(median = median(intensity)), by = get(compare_between)]$median)
  featureVals[, normIntensity := normIntensity * medianOfMedians]
  return(featureVals)
}





