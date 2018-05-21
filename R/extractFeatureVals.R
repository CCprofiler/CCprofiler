#' Extract intensity, correlation etc. within feature boundaries
#' @description Goes through protein/complex features and extracts metrics from the traces object
#' @param traces object of class traces
#' @param features Protein/complex feature table returned by \code{findProteinFeatures} or
#' \code{findComplexFeatures}.
#' @param perturb_cutoff Character string or numeric, the quantile of values that noise is sampled from.
#' Noise needs to be imputed to calculate the correlations.
#' @param verbose Logical, whether to print messages to console.
#' @param extract Character string, name ot the feature table column containing the subunits to extract.
#' @param imputeZero Logical, whether zero values of a fraction within a feature should be imputed, default = FALSE.
#' @return Long format table containing extracted feature values.
#' @export

extractFeatureVals <- function(traces, features,
                               perturb_cutoff = "5%",
                               verbose = TRUE,
                               extract = "subunits_detected",
                               imputeZero = FALSE,
                               fill = FALSE, ...){
  UseMethod("extractFeatureVals", traces)
}


#' @describeIn extractFeatureVals Extract values from single traces object.
#' @export
extractFeatureVals.traces <- function(traces, features,
                               perturb_cutoff = "5%",
                               verbose = TRUE,
                               extract = "subunits_detected",
                               imputeZero = FALSE, ...){

  if(!is.data.table(features)){
    stop("features must be of type 'data.table'")
  }
  featureColnames <- names(features)
  if("complex_id" %in% featureColnames){
    complexLevel <- TRUE
    if (traces$trace_type == "protein") {
      .tracesTest(traces, type = "protein")
      tracesType <- "protein"
      allIds <- features$complex_id
    } else if (traces$trace_type == "peptide") {
      tracesType <- "peptide"
      .tracesTest(traces, type = "peptide")
      allIds <- features$complex_id
    } else {
      stop("Could not detect if traces object is of type peptide or protein. Aborting...")
    }
  } else if("protein_id" %in% featureColnames){
    complexLevel <- FALSE
    if (traces$trace_type == "peptide") {
      tracesType <- "peptide"
      .tracesTest(traces, type = "peptide")
      allIds <- features$protein_id
    } else {
      stop("Could not detect if traces object is of type peptide. Aborting...")
    }
  } else {
    stop("Could not detect if feature table is for protein or complex Features. Aborting...")
  }

  #if(!("peak_corr" %in% featureColnames)) stop("No column peak_corr found. peak_corr must be
  #                                                calculated with calculateFeatureCorrelation() first!")
  if(!("peak_corr" %in% featureColnames)){features[,peak_corr:=1]} # I've done this to easily enable running feature extraction on collapsed complex features.

  # Get traces matrix
  traceMat <- getIntensityMatrix(traces)
  # Impute noise for missing intensity measurements globally for all traces
  traceMatImputed <- copy(traceMat)
  nZero <- sum(traceMatImputed == 0) # number of ZERO values in matrix
  measuredVals <- traceMatImputed[traceMatImputed != 0]
  if(class(perturb_cutoff) == "character"){
    qt <- as.numeric(gsub("%","",perturb_cutoff))/100
    perturb_cutoff <- quantile(measuredVals, qt)
  }
  set.seed(123) # set seed to always get same results
  traceMatImputed[traceMatImputed == 0] <- runif(min = 0, max = perturb_cutoff, nZero)

  # Calculate correlation in Peak boundaries for every detected trace

  featureVals <- apply(features, 1, function(hyp){

    # featuresHyp <- features[features$protein_id == hyp,]
    subunitsUnion <- unique(strsplit(hyp[extract],";")[[1]])
    if (complexLevel == FALSE) {
      subunitsUnion <- subunitsUnion[subunitsUnion %in% rownames(traceMat)]
    } else if ((complexLevel == TRUE) & (tracesType == "protein")) {
      subunitsUnion <- subunitsUnion[subunitsUnion %in% rownames(traceMat)]
    } else if ((complexLevel == TRUE) & (tracesType == "peptide")) {
      prot_ids <- traces$trace_annotation$protein_id
      subunitsUnion <- subunitsUnion[subunitsUnion %in% prot_ids]
      # change subunitsUnion to peptide ids
      subunitsUnion <- traces$trace_annotation[protein_id %in% subunitsUnion]$id
      pepProtMap <- subset(traces$trace_annotation[id %in% subunitsUnion],select=c("id","protein_id"))
      names(pepProtMap) <- c("pep_id","prot_id")
    } else {
      stop("Features and traces do not match. Aborting...")
    }
    nSubunits <- length(subunitsUnion)
    if(complexLevel){
      id <- as.character(hyp["complex_id"])
      name <- as.character(hyp["complex_name"])
    }else{
      id <- as.character(hyp["protein_id"])
      name <- as.character(hyp["protein_name"])
    }
    bound_left <- as.numeric(hyp["left_pp"])
    bound_right <- as.numeric(hyp["right_pp"])
    apex <- as.numeric(hyp["apex"])
    pk <- as.numeric(hyp["peak_corr"])
    # get intensity of every trace within peak boundaries
    tracesFeature <- traceMat[subunitsUnion,bound_left:bound_right]
    tracesFeatureImputed <- traceMatImputed[subunitsUnion,bound_left:bound_right]
    tracesAll <- traceMat[subunitsUnion,]
    tracesAllImputed <- traceMatImputed[subunitsUnion,]

    if(nSubunits>1){
        corr <- cor(t(tracesFeatureImputed))
        # get mean correlation for every subunit
        corrSubunits <- (colSums(corr)-1) / ((nSubunits-1)) #*pk to normalize?
        intensities <- rowSums(tracesFeature)
        rank <- rank(intensities)
        totalIntensity <- sum(intensities)
        totalTop2Intensity <- sum(sort(intensities, decreasing = T)[1:2])
        globalIntensity <- rowSums(tracesAll)
        globalIntensityImputed <- rowSums(tracesAllImputed)
        res <- data.table(tracesFeature)

    } else if(nSubunits == 1){
      if(verbose){
        message(paste0("Feature with only one subunit detected: ", id, " - Apex", apex,
                       ". This will produce NAs"))
      }
      corrSubunits = NA
      intensities <- sum(tracesFeature)
      rank <- 1
      totalIntensity <- intensities
      totalTop2Intensity <- NA
      globalIntensity <- sum(tracesAll)
      globalIntensityImputed <- sum(tracesAllImputed)
      res <- setDT(as.list(tracesFeature))[]
      # names(res) <- as.character(bound_left:bound_right)
    }else{
      if(verbose){
        message(paste0("Feature: ", id, " - Apex", apex, " has no detected subunits. Omitting."))
      }
      return(NULL)
    }

    if (tracesType=="peptide"){
      if (complexLevel){
        info_cols <- c("id","feature_id", "complex_id", "apex", "bound_left", "bound_right",
                       "corr","peak_cor", "total_pep_intensity", "total_complex_intensity",
                       "total_top2_complex_intensity", "global_intensity", "global_intensity_imputed")
        prot_id <- unlist(lapply(subunitsUnion,function(x){pepProtMap[pep_id==x]$prot_id}))
        res <- res[, (info_cols) := .(subunitsUnion, prot_id, id, apex, bound_left, bound_right, corrSubunits,
                                     pk, intensities, totalIntensity, totalTop2Intensity, globalIntensity, globalIntensityImputed)]
      } else {
        info_cols <- c("id","feature_id", "apex", "bound_left", "bound_right",
                       "corr","peak_cor", "total_pep_intensity", "total_prot_intensity",
                       "total_top2_prot_intensity", "global_intensity", "global_intensity_imputed")
        res <- res[, (info_cols) := .(subunitsUnion, id, apex, bound_left, bound_right, corrSubunits,
                                       pk, intensities, totalIntensity, totalTop2Intensity, globalIntensity, globalIntensityImputed)]
      }
    } else {
      info_cols <- c("id","feature_id", "apex", "bound_left", "bound_right",
                     "corr","peak_cor", "total_prot_intensity", "total_complex_intensity",
                     "total_top2_complex_intensity", "global_intensity", "global_intensity_imputed")
      res <- res[, (info_cols) := .(subunitsUnion, id, apex, bound_left, bound_right, corrSubunits,
                                   pk, intensities, totalIntensity, totalTop2Intensity, globalIntensity, globalIntensityImputed)]
    }
    resLong <- melt(res, id.vars = info_cols, variable.name = "fraction", value.name = "intensity")
    resLong
  })

  featureVals <- do.call("rbind", featureVals)
  featureVals <- as.data.table(featureVals)
  featureVals[, fraction := as.numeric(levels(fraction))[fraction]]

  ## Zero Value imputation
  if (imputeZero) {
    if(verbose) message("Filling in fractions with zero values within features...")
    set.seed(123)
    featureVals[, imputedFraction := (intensity == 0)]
    ## Completely absent features are imputed below the minimum value of that trace
    featureVals[, imputedFeature := (sum(intensity) == 0), by=.(id, apex)]
    if (TRUE %in% unique(featureVals$imputedFeature)) {
      featureVals[imputedFeature == TRUE,
                  min_int := unlist(lapply(id, function(x) min(traceMat[x,][traceMat[x,]>0])))]
      featureVals[imputedFeature == TRUE, intensity := runif(min = 0, max = min_int, .N)]
    }
    ## Features with some zeroes are imputed below the smallest value of that feature
    featureVals[intensity==0]$intensity <- NA
    featureVals[imputedFeature == FALSE, min_int := min(intensity,na.rm=TRUE), by=c("id","feature_id","apex")]
    ## featureVals[,new := unlist(lapply(min_int,function(x){sample(x,1)}))]
    ## featureVals[,intensity := ifelse(imputedFraction == TRUE, new, intensity)]
    featureVals[is.na(intensity), intensity := runif(min = 0, max = min_int, .N)]
    featureVals[, min_int := NULL]
    ## featureVals[, new := NULL]
  } else {
    featureVals[, imputedFraction := FALSE]
  }
  return(featureVals)
}

#' @describeIn extractFeatureVals Extract values from single traces object.
#' @export
extractFeatureVals.tracesList <- function(traces, features,
                                      perturb_cutoff = "5%",
                                      verbose = TRUE,
                                      design_matrix = NULL,
                                      extract = "subunits_detected",
                                      imputeZero = FALSE,
                                      fill = FALSE, ...){
  res <- lapply(names(traces), function(tr){
    message(paste0("Extracting values from ", tr))
    vals <- extractFeatureVals.traces(traces[[tr]],
                                      features = features,
                                      perturb_cutoff = perturb_cutoff,
                                      verbose = verbose,
                                      extract = extract,
                                      imputeZero = imputeZero)
    vals[,Sample := tr]
    return(vals)
  })
  res <- do.call(rbind, res)
  ## Remove any traces that are found in no condition
  res <- res[res[, .I[!all(imputedFraction)], by=c("id","feature_id","apex")]$V1]

  if(!is.null(design_matrix)){
    res <- merge(res, unique(design_matrix[Sample_name %in% unique(res$Sample), .(Sample_name, Condition)]),
                 by.x = "Sample", by.y = "Sample_name")
    res <- merge(res, unique(design_matrix[Sample_name %in% unique(res$Sample), .(Sample_name, Replicate)]),
                 by.x = "Sample", by.y = "Sample_name")

    if(fill){
      if(verbose) message("Filling in missing values across Features...")
      res <- fillFeatureVals(featureVals = res,
                                     design_matrix = design_matrix,
                                     perturb_cutoff = perturb_cutoff)
    }
  }else if(fill){
    message("Cannot fill values without a design matrix. Please specify.")
    message("Will return non-filled feature Values...")
  }
  return(res)
}

#' Fill traceVals for all Samples
#' @description Fill in traces that are not detected in certain conditions with noise values
#' @param featureVals data.table, a long-format table of intensity values within feature boundaries.
#' A featureVals table can be produced with \code{extractFeatureVals}.
#' @param design_matrix data.table, design matrix describing the architecture of the tracesList object.
#' @param perturb_cutoff Character string or numeric, the quantile of values that noise is sampled from.
#' Noise needs to be imputed to calculate the correlations.
#' @return Long format table containing filled feature values.
#' @export

fillFeatureVals <- function(featureVals,
                            design_matrix,
                            perturb_cutoff = "5%"){

  if(class(perturb_cutoff) == "character"){
    qt <- as.numeric(gsub("%","",perturb_cutoff))/100
    perturb_cutoff <- quantile(featureVals$intensity, qt)
  }
  set.seed(123) # set seed to always get same results
  samples <- unique(design_matrix$Sample)
  n_samples <- length(samples)
  fv <- copy(featureVals)
  if ("complex_id" %in% names(fv)) {
    key_v <- c("id", "feature_id", "complex_id", "apex", "fraction", "bound_left", "bound_right", "Sample")
    setkeyv(fv, key_v)
    complete_table <- unique(fv[, .(id, feature_id, complex_id, apex, fraction, bound_left, bound_right)])
  } else {
    key_v <- c("id", "feature_id", "apex", "fraction", "bound_left", "bound_right", "Sample")
    setkeyv(fv, key_v)
    complete_table <- unique(fv[, .(id, feature_id, apex, fraction, bound_left, bound_right)])
  }

  # complete_table_ <- apply(complete_table, 1, function(x) cbind(rbind(rep(x, n_samples)), samples))
  complete_table_ <- do.call(rbind,lapply(samples, function(x) cbind(complete_table, x)))
  setnames(complete_table_, "x", "Sample")
  setkeyv(complete_table_, key(fv))
  fvComp <- merge(complete_table_, fv,
                  by=key_v,
                  all.x= T)
  fvComp[is.na(imputedFraction)]$imputedFraction <- TRUE
  fvComp[is.na(imputedFeature)]$imputedFeature <- TRUE
  ## Remove any traces that are found in no condition
  fvComp <- fvComp[fvComp[, .I[!all(imputedFraction)], by=c("id","feature_id","apex")]$V1]

  fvComp[, imputedCondition := is.na(intensity)]
  n_filled <- fvComp[imputedCondition == T, .N]
  ## Impute missing features with the minimum Value with the perturbation cutoff
  fvComp[imputedCondition == T, intensity := as.double(runif(min = 0, max = perturb_cutoff, n_filled))]
  ## Impute missing traces with the feature minimum of the other fraction
  fvComp[, min_int := min(intensity[imputedFraction == F & intensity > 0],na.rm=TRUE), by=c("id","feature_id","apex")]
  fvComp[imputedCondition  == T, intensity := runif(min=0, max = min_int, n=n_filled)]

  fvComp[, Condition := NULL]
  fvComp[, Replicate := NULL]
  fvComp <- merge(fvComp, design_matrix[,.(Sample_name, Condition)],
                  by.x = "Sample", by.y = "Sample_name", sort = F)
  fvComp <- merge(fvComp, design_matrix[,.(Sample_name, Replicate)],
                  by.x = "Sample", by.y = "Sample_name", sort = F)

  return(fvComp[])
}
