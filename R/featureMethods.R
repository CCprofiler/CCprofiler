#' Summarize co-elution features
#' @description Summarize co-elution features
#' @param feature_table data.table as reported by \code{\link[SECprofiler]{findComplexFeatures}} 
#' or \code{\link[SECprofiler]{findProteinFeatures}}.
#' @param plot Logical, whether to generate a feature summary plot. Default is \code{TRUE}.
#' @param PDF Logical, whether to generate a PDF file with the summary plot. Default is \code{FALSE}.
#' @param name Character string specifying the name of the PDF file of the summary plot.
#' Only applicable if PDF=\code{TRUE}. Default is "feature_summary".
#' @return A summary list with complex feature metrics.
#' @examples
#' ## Load example complex feature finding results:
#' complexFeatures <- exampleComplexFeatures
#' ## Run summary function:
#' summarizeFeatures(complexFeatures)
#' ## Filter complex features:
#' filteredComplexFeatures <- filterFeatures(complexFeatures)
#' ## Run summary function on filtered data:
#' summarizeFeatures(filteredComplexFeatures)
#' @export
summarizeFeatures <- function(feature_table,
                              plot=TRUE,
                              PDF=FALSE,
                              name="feature_summary"){
  features <- copy(feature_table)
  if ("protein_id" %in% names(features)){
    type <- "protein"
    setnames(features,"protein_id","complex_id")
  } else if ("complex_id" %in% names(features)){
    type <- "complex"
  } else {
    stop("This is no protein or complex feature table. Please provide features from findComplexFeatures or findProteinFeatures.")
  }
  totalConfirmedHypotheses <- length(unique(features$complex_id))
  totalFeatures <- nrow(features)
  features[ , `:=`( COUNT = .N , IDX = 1:.N ) , by = complex_id ]
  # setkey(features, "complex_id")
  feature_count_max <- unique(features, by = "complex_id")$COUNT
  summaryFeatureCount <- summary(feature_count_max) 
  totalHypothesesWithMultipleFeatures <- length(which(feature_count_max>=2))
  summaryCorrelation <- summary(features$peak_corr)
  summaryArea <- summary(features$area)
  if("mw_diff" %in% names(features)) {
    summaryMWdiff <- summary(features$mw_diff)
  }
  summaryNsubunitsAnnotated <- summary(features$n_subunits_annotated)
  summaryNsubunitsWithSignal <- summary(features$n_subunits_with_signal)
  summaryNsubunitsDetected <- summary(features$n_subunits_detected)
  summaryCompleteness <- summary(features$completeness)
  if("mw_diff" %in% names(features)) {
    summary <- list(type = type,
                   totalFeatures = totalFeatures,
                   totalConfirmedHypotheses = totalConfirmedHypotheses,
                   totalHypothesesWithMultipleFeatures = totalHypothesesWithMultipleFeatures,
                   summaryFeatureCount = summaryFeatureCount,
                   summaryCorrelation = summaryCorrelation,
                   summaryArea = summaryArea,
                   summaryMWdiff = summaryMWdiff,
                   summaryNsubunitsAnnotated = summaryNsubunitsAnnotated,
                   summaryNsubunitsWithSignal = summaryNsubunitsWithSignal,
                   summaryNsubunitsDetected = summaryNsubunitsDetected,
                   summaryCompleteness = summaryCompleteness
                   )
  } else {
    summary <- list(type = type,
                    totalFeatures = totalFeatures,
                    totalConfirmedHypotheses = totalConfirmedHypotheses,
                    totalHypothesesWithMultipleFeatures = totalHypothesesWithMultipleFeatures,
                    summaryFeatureCount = summaryFeatureCount,
                    summaryCorrelation = summaryCorrelation,
                    summaryArea = summaryArea,
                    summaryNsubunitsAnnotated = summaryNsubunitsAnnotated,
                    summaryNsubunitsWithSignal = summaryNsubunitsWithSignal,
                    summaryNsubunitsDetected = summaryNsubunitsDetected,
                    summaryCompleteness = summaryCompleteness
    )
  }
  if (plot) {
    features[,feature_id := .I]
    if("mw_diff" %in% names(features)) {
      plottingFeatures <- subset(features,select=c("feature_id","peak_corr","mw_diff","area","apex","n_subunits_detected","completeness")) 
    } else {
      plottingFeatures <- subset(features,select=c("feature_id","peak_corr","area","apex","n_subunits_detected","completeness")) 
    }
    plottingFeatures$area <- log(plottingFeatures$area)
    setnames(plottingFeatures,"area","log(area)")
    plottingFeatures[, names(plottingFeatures) := lapply(.SD, as.numeric)]
    plottingFeaturesMelt <- data.table::melt(plottingFeatures,id.vars=c("feature_id"))
    if (PDF) {
      pdf(gsub("$|\\.pdf$", ".pdf", name))
    }
    p <- ggplot(plottingFeaturesMelt, aes(x=value,fill=variable)) + 
      geom_histogram(bins = 50) + 
      facet_wrap(~variable, scales="free", ncol = 2) + 
      theme_bw() +
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + 
      guides(fill=FALSE) +
      ggtitle(paste0(type," feature summary")) + 
      theme(plot.title = element_text(size=14, face="bold"))
    print(p)
    plottingNsubcomplexes <- data.table(n_subcomplexes=feature_count_max)
    q <- ggplot(plottingNsubcomplexes,aes(x=n_subcomplexes)) +
      stat_bin(binwidth=1) +
      stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-0.5) +
      labs(x="N co-elution features",y="N complex hypotheses") +
      theme_classic() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      ggtitle(paste0(type," sub-feature summary")) + 
      theme(plot.title = element_text(size=14, face="bold"))
    print(q)
    if (PDF) {
      dev.off()
    }
  }
  return(summary)
}

#' Select best feature for each complex or protein id.
#' @description Filter feature table to only contain one best feature per complex or protein id.
#' "Best" is defined as the feature with most subunits, highest peak correlation and largest signal area.
#' @param feature_table data.table as reported by \code{\link[SECprofiler]{findComplexFeatures}} 
#' or \code{\link[SECprofiler]{findProteinFeatures}}.
#' @return data.table Filtered feature table with only one best feature per complex or protein id.
#' @examples
#' ## Load example complex feature finding results:
#' complexFeatures <- exampleComplexFeatures
#' ## Run summary function:
#' summarizeFeatures(complexFeatures)
#' ## Filter complex features to only contain one best feature per complex id:
#' bestComplexFeatures <- getBestFeatures(complexFeatures)
#' ## Run summary function on filtered data:
#' summarizeFeatures(bestComplexFeatures)
#' @export
getBestFeatures <- function(feature_table){
  if ("complex_id" %in% names(feature_table)) {
    feature_table <- feature_table[order(complex_id,-n_subunits_detected,-peak_corr,-area)]
    res_best <- unique(feature_table,by="complex_id")
    return(res_best)
  } else if ("protein_id" %in% names(feature_table)) {
    feature_table <- feature_table[order(protein_id,-n_subunits_detected,-peak_corr,-area)]
    res_best <- unique(feature_table,by="protein_id")
    return(res_best)
  } else {
    stop("This is no protein or complex feature table. Please provide features from findComplexFeatures or findProteinFeatures.")
  }
}


#' Filter co-elution feature table
#' @description Filter co-elution feature table according to desired criteria.
#' @param feature_table data.table as reported by \code{\link[SECprofiler]{findComplexFeatures}} 
#' or \code{\link[SECprofiler]{findProteinFeatures}}.
#' @param complex_ids character vector containing all desired \code{complex_id} values. 
#' Only applicable to complex feature tables. Defaults to \code{NULL}.
#' @param protein_ids character vector containing all desired \code{protein_id} values. 
#' Only applicable to protein feature tables. Defaults to \code{NULL}.
#' @param min_feature_completeness Numeric between 0 and 1, specifying the required completeness 
#' of a feature reltive to the tested hypothesis (keeps all features if at least one is bigger than the cutoff).
#' @param min_hypothesis_completeness Numeric between 0 and 1, specifying the required completeness 
#' of the most complete feature reltive to the tested hypothesis (keeps all features if at least one feature is bigger than the cutoff).
#' @param min_subunits Integer specifying minimum number of subunits in a co-elution feature.
#' @param min_peak_corr Numeric value betwee 0 and 1 specifying minimum peak correlation, defaults to 0.5.
#' @param min_monomer_distance_factor Numeric value specifying a factor to multiply the largest monomer molecular weight, defaults to \code{NULL}.
#' This filters out features that have their apex at a smaller molecular weight than the resulting min_monomer_distance_factor*max(monomer_mw) value.
#' Using this filtering option can, for example, remove complex features that are likely spontaneous co-elutions of the subunits monomers.
#' @return The same feture table as teh input, but filtered according to the provided parameters.
#' @examples
#' ## Load example complex feature finding results:
#' complexFeatures <- exampleComplexFeatures
#' ## Run summary function:
#' summarizeFeatures(complexFeatures)
#' ## Filter complex features by a peak correlation of 0.5, a minimum hypothesis
#' ## completeness of 0.5 and a minimum distance to the monomers by a factor of 2:
#' filteredComplexFeatures <- filterFeatures(complexFeatures, 
#'                                      min_peak_corr=0.5, 
#'                                      min_hypothesis_completeness=0.5,
#'                                      min_monomer_distance_factor=2)
#' ## Run summary function on filtered data:
#' summarizeFeatures(filteredComplexFeatures)
#' @export
filterFeatures <- function(feature_table,
                           complex_ids=NULL,
                           protein_ids=NULL,
                           min_feature_completeness=NULL,
                           min_hypothesis_completeness=NULL,
                           min_subunits=NULL,
                           min_peak_corr=0.5,
                           min_monomer_distance_factor=NULL){
  if("complex_id" %in% names(feature_table)){
    type="complex"
  } else if ("protein_id" %in% names(feature_table)) {
    type="protein"
  } else {
    stop("This is neither a complex or a protein feature table. Please check the input data.")
  }
  if(!is.null(complex_ids)){
    if(type=="complex"){
      feature_table <- subset(feature_table,complex_id %in% as.character(complex_ids))
    } else {
      stop("This is not a complex feature table. Please check the input data and your filtering criteria, especially the complex_ids.")
    }
  }
  if(!is.null(protein_ids)){
    if(type=="protein"){
      feature_table <- subset(feature_table,protein_id %in% as.character(protein_ids))
    } else {
      stop("This is not a protein feature table. Please check the input data and your filtering criteria, especially the protein_ids.")
    }
  }
  if(!is.null(min_feature_completeness)){
    feature_table <- subset(feature_table,completeness >= min_feature_completeness)
  }
  if(!is.null(min_hypothesis_completeness)){
    if(type=="complex"){
      allowed_ids <- feature_table[completeness>=min_hypothesis_completeness, unique(complex_id)]
      feature_table <- subset(feature_table,complex_id %in% allowed_ids)
    } else if (type=="protein") {
      allowed_ids <- feature_table[completeness>=min_hypothesis_completeness, unique(protein_id)]
      feature_table <- subset(feature_table,protein_id %in% allowed_ids)
    }
  }
  if(!is.null(min_subunits)){
    feature_table <- subset(feature_table,n_subunits_detected >= min_subunits)
  }
  if(!is.null(min_peak_corr)){
    feature_table <- subset(feature_table,peak_corr >= min_peak_corr)
  }
  if(!is.null(min_monomer_distance_factor)){
    if (! "monomer_mw" %in% names(feature_table)){
      stop("No molecular weight information available. The min_monomer_distance_factor is only a valid 
           filtering option if molecular weight information of the proteins is available.")
    }
    dist <- lapply(seq(1:nrow(feature_table)), function(i){
      feature=feature_table[i]
      max_monomer_mw <- max(as.numeric(strsplit(as.character(feature$monomer_mw), ';')[[1]]))
      mw_dist <- feature$apex_mw-(min_monomer_distance_factor*max_monomer_mw)
      mw_dist
    })
    dist_vector <- unlist(dist)
    sel <- which(dist<=0)
    if(length(sel)>0){
      feature_table <- feature_table[-sel]
    }
  }
  return(feature_table)
}


#' Filter co-elution feature table by stepwise completeness cutoffs
#' @description Filter co-elution feature table by two consecutive completeness cutoffs 
#' to treat small and large complexes differenetly.
#' @param feature_table data.table as reported by \code{\link[SECprofiler]{findComplexFeatures}} 
#' or \code{\link[SECprofiler]{findProteinFeatures}}.
#' @param min_subunits_annotated Integer specifying the number of annotated hypothesis components (peptides / proteins).
#' This is the cutoff for using the two different completeness cutoffs as provided by the \code{completeness_vector}. Default is \code{4}.
#' @param completeness_vector Numeric vector of length 2. The first value is the completeness cutoff apllied to all hypotheses
#' with <= \code{min_subunits_annotated} subunits and the secind value is applied to all hypotheses 
#' with more than \code{min_subunits_annotated} subunits. Default is \code{c(0.75,0.5)}
#' @param level Character string defining level of filterng, allowed values are "feature" or "hypothesis". 
#' "feature" filters all features by their completeness. "hypothesis" filters by the completeness of the 
#' of the largest feature within one hypothesis. Default is "hypothesis".
#' @return The same feture table as teh input, but filtered according to the provided parameters.
#' @examples
#' ## Load example complex feature finding results:
#' complexFeatures <- exampleComplexFeatures
#' ## Run summary function:
#' summarizeFeatures(complexFeatures)
#' ## Filter complex features by a subunit cutoff of 4 and the completness cutoffs
#' ## 0.75 for hypotheses <= 4 subunits and 0.5 for hypotheses with > 4 annotated subunits.
#' filteredComplexFeatures <- filterByStepwiseCompleteness(feature_table=complexFeatures,
#'                                           min_subunits_annotated=4,
#'                                           completeness_vector=c(0.75,0.5),
#'                                           level="hypothesis")
#' ## Run summary function on filtered data:
#' summarizeFeatures(filteredComplexFeatures)
#' @export
filterByStepwiseCompleteness <- function(feature_table,
                                         min_subunits_annotated=4,
                                         completeness_vector=c(0.75,0.5),
                                         level="hypothesis"){
  if (length(completeness_vector) != 2) {
    stop("The completeness vector should contain m=2 values.")
  }
  features <- copy(feature_table)
  if ("protein_id" %in% names(features)){
    type <- "protein"
    setnames(features,"protein_id","complex_id")
  } else if ("complex_id" %in% names(features)){
    type <- "complex"
  } else {
    stop("This is no protein or complex feature table. Please provide features from findComplexFeatures or findProteinFeatures.")
  }
  subset_small <- subset(features, n_subunits_annotated <= min_subunits_annotated)
  subset_big <- subset(features, n_subunits_annotated > min_subunits_annotated)
  if (level=="hypothesis") {
    allowed_ids_small <- subset_small[completeness>=completeness_vector[1],unique(complex_id)]
    allowed_ids_big <- subset_big[completeness>=completeness_vector[2],unique(complex_id)]
    features_filtered <- subset(features,complex_id %in% c(allowed_ids_small,allowed_ids_big))
    return(features_filtered)
  } else if (level=="feature") {
    subset_small <- subset(subset_small,completeness>=completeness_vector[1])
    subset_big <- subset(subset_big,completeness>=completeness_vector[2])
    subset <- rbind(subset_small,subset_big)
    return(subset)
  } else {
    stop("No valid level provided. The oly valid levels are \"hypothesis\" or \"feature\".")
  }
}


