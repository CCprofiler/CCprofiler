#' A sample data.table that holds information on CORUM complex membership.
#'
#' Each row in this data.table corresponds to an association between a CORUM
#' complex and a protein subunit.
#' The structure of this data.table is as follows:
#' \itemize{
#'  \item \code{complex_id} The numeric id of a CORUM complex.
#'  \item \code{complex_name} A character string of the complex' name.
#'  \item \code{protein_id} The Uniprot identifier of the protein.
#' }
'corumComplexHypotheses'

#' A sample data.table that holds information on example protein complex hypotheses.
#'
#' Each row in this data.table corresponds to an association between a CORUM
#' complex and a protein subunit.
#' The structure of this data.table is as follows:
#' \itemize{
#'  \item \code{complex_id} The numeric id of a complex.
#'  \item \code{complex_name} A character string of the complex' name.
#'  \item \code{protein_id} The Uniprot identifier of the protein.
#' }
'exampleComplexHypotheses'

#' A sample data.table representing the output format of the OpenSWATH workflow.
'exampleOpenSWATHinput'

#' A sample data.table that holds information on example protein complex hypotheses.
#'
#' Each row in this data.table corresponds to an association between a filename
#' and the SEC fraction. The fraction_numbers should span a range from 1 to
#' the total number of fractions, without duplicates.
#' The structure of this data.table is as follows:
#' \itemize{
#'  \item \code{filename} Filename as in the input proteomic data.
#'  \item \code{fraction_number} Fraction number assotiated to the filename.
#' }
'exampleFractionAnnotation'

#' A sample data.table that holds information on proteins or peptides.
#'
#' Each row in this data.table corresponds to an association between a protein/peptide
#' and any specific information, for example annotations from Uniprot.
'exampleTraceAnnotation'

#' A sample peptide level traces object before any filtering.
'examplePeptideTraces'

#' A sample peptide level traces object after consecutive id and sibling peptide correlation filter.
'examplePeptideTracesFiltered'

#' A sample protein level traces object.
'exampleProteinTraces'

#' A sample data.table that holds information on the connectivity between SEC fraction number and apparent
#' molecular weight (in kDa).
#'
#' Each row in this data.table corresponds to an association between SEC fraction number and apparent
#' molecular weight (in kDa).
#' The structure of this data.table is as follows:
#' \itemize{
#'  \item \code{std_weights_kDa} The molecular weights of the standard proteins.
#'  \item \code{std_elu_fractions} The fraction numbers where these standard proteins eluted.
#' }
'exampleCalibrationTable'
