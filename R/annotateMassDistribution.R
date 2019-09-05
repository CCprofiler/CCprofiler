#' annotateMassDistribution
#' @description Estimates fraction of assembled vs. monomeric mass per id.
#' @param traces An object of type traces.
#' @export
annotateMassDistribution <- function(traces){
  .tracesTest(traces)

  if(! ("protein_mw" %in% names(traces$trace_annotation))) {
    stop("annotateMassDistribution requires information on the monomer molecular weight. Please see annotateTraces.")
  }

  if(! ("molecular_weight" %in% names(traces$fraction_annotation))) {
    stop("annotateMassDistribution requires information on the molecular weight assigned to each fraction. Please see annotateMolecularWeight.")
  }

  intMat <- getIntensityMatrix(traces)
  mw_info <- traces$trace_annotation
  max_trace <- max(traces$fraction_annotation$id)
  mw_info[,assembly_boundary := 2*protein_mw]
  mw_info[is.na(assembly_boundary), assembly_boundary := 0]
  mw_info[,assembly_boundary_fraction := unlist(lapply(assembly_boundary,function(x){max(traces$fraction_annotation[molecular_weight>=x]$id)}))]

  # Get intensity sum
  mw_info[,sum_assembled := sum(intMat[id,1:assembly_boundary_fraction-1]), id]
  mw_info[,sum_monomeric := sum(intMat[id,assembly_boundary_fraction:max_trace]), id]
  mw_info[,sum := sum(intMat[id,]), id]
  mw_info[is.na(assembly_boundary), sum_assembled := 0]
  mw_info[is.na(assembly_boundary), sum_monomeric := 0]
  mw_info[, sum_assembled_norm := 1/sum*sum_assembled]
  mw_info[, sum_monomeric_norm := 1/sum*sum_monomeric]
  #mw_info[, test := sum_monomeric_norm+sum_assembled_norm]

  traces$trace_annotation <- mw_info
  .tracesTest(traces)
  return(traces)
}

#' summarizeMassDistribution
#' @description Plots the fraction of assembled vs. monomeric mass
#' @param traces An object of type traces.
#' @param PDF Logical, whether to generate a PDF file with the summary plot. Default is \code{FALSE}.
#' @param name Character string specifying the name of the PDF file of the summary plot.
#' Only applicable if PDF=\code{TRUE}. Default is "massDistribution_summary".
#' @export
summarizeMassDistribution <- function(traces,
                              PDF=FALSE,
                              name="massDistribution_summary"){
  .tracesTest(traces)
  if(traces$trace_type != "protein") {
    stop("Traces need to be of type protein for this function.")
  }
  protTraces_annotated <- annotateMassDistribution(traces)
  monomer_mass <- sum(protTraces_annotated$trace_annotation$sum_monomeric)
  assembled_mass <- sum(protTraces_annotated$trace_annotation$sum_assembled)
  full_mass <- sum(protTraces_annotated$trace_annotation$sum)

  if (PDF) {
    pdf(gsub("$|\\.pdf$", ".pdf", name))
   }
 pie_plot <- pie(c(monomer_mass,assembled_mass),
                labels=c(paste0("in monomer range ",round(monomer_mass/full_mass, 2) *100,"%"),
                         paste0("in assembled range ",round(assembled_mass/full_mass, 2) *100,"%")),
                init.angle = 180, main = "Total MS signal \n(~ protein mass)")
  print(pie_plot)
  if (PDF) {
    dev.off()
   }
}
