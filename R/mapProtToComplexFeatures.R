#' Map protein features to complex features
#' @description Map co-elution features from the protein centric feature finding
#' to co-elution features from the complex-centric workflow.
#' @import data.table
#' @param proteinFeatures data.table with protein co-elution features in the
#' format of the output from \code{proteinFeatureFinding} or
#' \code{aggregatePeptideTests}.
#' @param complexFeatures data.table with complex co-elution features in the
#' format of the output from \code{complexFeatureFinding}.
#' @param prot_id chracter string specifying the column name of the id type
#' in the proteinFeatures data.table to match, default is "protein_id".
#' @param max_apex_dist numeric specifying the minimum number of fractions
#' between a protein and complex co-elution feature to be mapped 2, default is 2.
#' @param plot Logical if histogram of apex distances between protein and
#' complex co-elution features should be plotted, default is FALSE.
#' @param PDF Logical if plot should be written to PDF, default is FALSE.
#' @return data.table with protein co-elution features that have a matching
#' complex co-elution feature.
#' @examples
#' # Load some example data:
#' proteinFeatures <- exampleProteinFeatures
#' complexFeatures <- exampleComplexFeatures
#' # Run the mapping:
#' mappedFeatures <- mapProtToComplexFeatures(proteinFeatures,
#'                                            complexFeatures,
#'                                            prot_id = "protein_id",
#'                                            max_apex_dist = 2
#'                                            )
#' # Inspect mapping result:
#' mappedFeatures
#' @export

mapProtToComplexFeatures <- function(proteinFeatures, complexFeatures, prot_id = "protein_id", max_apex_dist = 2, plot=F, PDF=F) {
  if (! all(c(prot_id,"apex") %in% names(proteinFeatures))) {
    stop(paste0("proteinFeatures do not contain minimum columns required: ",prot_id," and apex"))
  }
  if (! all(c("complex_id","subunits_detected","apex") %in% names(complexFeatures))) {
    stop("complexFeatures do not contain minimum columns required: complex_id, subunits_detected and apex")
  }
  complex_info <- subset(complexFeatures, select=c("complex_id","complex_name","subunits_detected","apex"))
  setnames(complex_info,c("complex_id","complex_name","subunits_detected","apex"),c("complex_id","complex_name","complex_subunits_detected","complex_apex"))
  complex_info <- complex_info[,list(protein_id = unlist(strsplit(complex_subunits_detected, ";"))), by=c("complex_id","complex_name","complex_subunits_detected","complex_apex")]
  mappedProtFeatures <- merge(proteinFeatures,complex_info,by.x=prot_id,by.y="protein_id",all.x=F,all.y=F)
  mappedProtFeatures[,apex_dist := abs(apex-complex_apex)]
  if (plot) {
    if (PDF) {
      pdf("mapping_dist_hist.pdf")
    }
      hist(mappedProtFeatures$apex_dist,breaks=seq(0,max(mappedProtFeatures$apex_dist),1))
    if (PDF) {
      dev.off()
    }
  }
  mappedProtFeatures <- mappedProtFeatures[apex_dist <= max_apex_dist]
  return(mappedProtFeatures)
}
