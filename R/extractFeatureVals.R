#' Extract intensity, correlation etc. within feature boundaries
#' @description Goes through protein/complex features and extracts metrics from the traces object
#' @param traces object of class traces
#' @param features Protein/complex feature table returned by \code{findProteinFeatures} or
#' \code{findComplexFeatures}.
#' @param perturb_cutoff Character string or numeric, the quantile of values that noise is sampled from.
#' Noise needs to be imputed to calculate the correlations.
#' @param verbose Logical, whether to print messages to console.
#' @param extract Character string, name ot the feature table column containing the subunits to extract.
#' @return Long format table containing extracted feature values.
#' @export

extractFeatureVals <- function(traces, features,
                               perturb_cutoff = "5%",
                               verbose = TRUE,
                               extract = "subunits_detected", ...){
  UseMethod("extractFeatureVals", traces)
}


#' @describeIn extractFeatureVals Extract values from single traces object.
#' @export
extractFeatureVals.traces <- function(traces, features,
                               perturb_cutoff = "5%",
                               verbose = TRUE, 
                               extract = "subunits_detected", ...){
  
  if(!is.data.table(features)){
    stop("features must be of type 'data.table'")
  }
  featureColnames <- names(features)
  if("complex_id" %in% featureColnames){
    complexLevel <- TRUE
    .tracesTest(traces, type = "protein")
    allIds <- features$complex_id
  } else if("protein_id" %in% featureColnames){
    complexLevel <- FALSE
    .tracesTest(traces, type = "peptide")
    allIds <- features$protein_id
  } else {
    stop("Could not detect if feature table is for protein or complex Features. Aborting...")
  }
  
  if(!("peak_corr" %in% featureColnames)) stop("No column peak_corr found. peak_corr must be
                                                  calculated with calculateFeatureCorrelation() first!")
  
  
  
  # Impute noise for missing intensity measurements globally for all traces
  traceMatImputed <- getIntensityMatrix(traces)
  nZero <- sum(traceMatImputed == 0) # number of ZERO values in matrix
  measuredVals <- traceMatImputed[traceMatImputed != 0]
  if(class(perturb_cutoff) == "character"){
    qt <- as.numeric(gsub("%","",perturb_cutoff))/100
    perturb_cutoff <- quantile(measuredVals, qt)
  }
  set.seed(123) # set seed to always get same results
  traceMatImputed[traceMatImputed == 0] <- sample(1:perturb_cutoff,size = nZero,
                                                      replace = TRUE)
  # Calculate correlation in Peak boundaries for every detected trace
  
  featureVals <- apply(features, 1, function(hyp){
    
    # featuresHyp <- features[features$protein_id == hyp,]
    subunitsUnion <- unique(strsplit(hyp[extract],";")[[1]])
    subunitsUnion <- subunitsUnion[subunitsUnion %in% rownames(traceMatImputed)]
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
    tracesFeature <- traceMatImputed[subunitsUnion,bound_left:bound_right]
    
    if(nSubunits>1){
        corr <- cor(t(tracesFeature))
        # get mean correlation for every subunit
        corrSubunits <- (colSums(corr)-1) / ((nSubunits-1)) #*pk to normalize?
        intensities <- rowSums(tracesFeature)
        rank <- rank(intensities)
        totalIntensity <- sum(intensities)
        totalTop2Intensity <- sum(sort(intensities, decreasing = T)[1:2])
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
      res <- setDT(as.list(tracesFeature))[]
      # names(res) <- as.character(bound_left:bound_right)
    }else{
      if(verbose){
        message(paste0("Feature: ", id, " - Apex", apex, " has no detected subunits. Omitting."))
      }
      return(NULL)
    } 
    
    info_cols <- c("id","feature_id", "apex", "bound_left", "bound_right",
                   "corr","peak_cor", "total_pep_intensity", "total_prot_intensity",
                   "total_top2_prot_intensity")
    res[, (info_cols) := .(subunitsUnion, id, apex, bound_left, bound_right, corrSubunits,
                           pk, intensities, totalIntensity, totalTop2Intensity)]
    resLong <- melt(res, id.vars = info_cols, variable.name = "fraction", value.name = "intensity")
  })
  featureVals <- do.call("rbind", featureVals)
  featureVals[, fraction := as.numeric(levels(fraction))[fraction]]
  return(featureVals)
}

#' @describeIn extractFeatureVals Extract values from single traces object.
#' @export
extractFeatureVals.tracesList <- function(traces, features,
                                      perturb_cutoff = "5%",
                                      verbose = TRUE,
                                      design_matrix = NULL,
                                      extract = "subunits_detected", ...){
  res <- lapply(names(traces), function(tr){
    message(paste0("Extracting values from ", tr))
    vals <- extractFeatureVals.traces(traces[[tr]], features, perturb_cutoff, verbose, extract)
    if(!is.null(design_matrix)){
      vals[,Condition := unique(design_matrix[Sample_name == tr, Condition])]
      vals[,Replicate := unique(design_matrix[Sample_name == tr, Replicate])]
    }
    vals[,Sample := tr]
    return(vals)
  })
  do.call(rbind, res)
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
  setkeyv(fv, c("id", "feature_id", "apex", "fraction", "bound_left", "bound_right", "Sample"))
  complete_table <- unique(fv[, .(id, feature_id, apex, fraction, bound_left, bound_right)])
  # complete_table_ <- apply(complete_table, 1, function(x) cbind(rbind(rep(x, n_samples)), samples))
  complete_table_ <- do.call(rbind,lapply(samples, function(x) cbind(complete_table, x)))
  setnames(complete_table_, "x", "Sample")
  setkeyv(complete_table_, key(fv))
  fvComp <- merge(complete_table_, fv,
                  by=c("id", "feature_id", "apex", "fraction", "bound_left", "bound_right", "Sample"),
                  all.x= T)
  fvComp[, filled_in := is.na(intensity)]
  n_filled <- fvComp[filled_in == T, .N]
  fvComp[filled_in == T, intensity := as.double(sample(1:perturb_cutoff, n_filled, replace = T))]
  return(fvComp)
}
  
  