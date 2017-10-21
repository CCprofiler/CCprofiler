#' Complex hypotheses from CORUM DB
#'
#' Each row in this data.table corresponds to an association between a CORUM
#' complex id/name and a protein subunit.
#' The structure of this data.table is as follows:
#' \itemize{
#'  \item \code{complex_id} The numeric id of a CORUM complex.
#'  \item \code{complex_name} A character string of the complex name.
#'  \item \code{protein_id} The Uniprot identifier of the protein.
#' }
'corumComplexHypotheses'

#' Complex hypotheses from CORUM DB (redundant)
#'
#' Each row in this data.table corresponds to an association between a CORUM
#' complex id/name and a protein subunit.
#' The structure of this data.table is as follows:
#' \itemize{
#'  \item \code{complex_id} The numeric id of a CORUM complex.
#'  \item \code{complex_name} A character string of the complex name.
#'  \item \code{protein_id} The Uniprot identifier of the protein.
#' }
'corumComplexHypothesesRedundant'

#' Example complex hypotheses
#'
#' Each row in this data.table corresponds to an association between a CORUM
#' complex id/name and a protein subunit. This is a reduced set of the 
#' \code{corumComplexHypotheses} data.table used for examples and testing. 
#' The structure of this data.table is as follows:
#' \itemize{
#'  \item \code{complex_id} The numeric id of a CORUM complex.
#'  \item \code{complex_name} A character string of the complex name.
#'  \item \code{protein_id} The Uniprot identifier of the protein.
#' }
'exampleComplexHypotheses'

#' Examplary annotation table that maps filenames to chromatographic fractions
#'
#' Each row in this data.table corresponds to an association between a filename
#' and a chromatographic fraction. The fraction_numbers should span a range from 1 to
#' the total number of fractions, without duplicates.
#' The structure of this data.table is as follows:
#' \itemize{
#'  \item \code{filename} Character string with the filename as in the input proteomic data.
#'  \item \code{fraction_number} Integer fraction number assotiated to the filename.
#' }
'exampleFractionAnnotation'

#' Examplary trace annotation table mapping ids in the proteomic data to Uniprot annotations
#'
#' Each row in this data.table corresponds to an association between a protein/peptide
#' and any specific information, for example annotations downloaded from Uniprot.
'exampleTraceAnnotation'

#' Examplary calibration table that maps molecular weights to elution fractions
#'
#' Each row in this data.table corresponds to an association between a molecular weight (kDa)
#' and a chromatographic fraction. These associations can be generated from standard proteins 
#' that are separated with the applied chromatographic separation technique.
#' The structure of this data.table is as follows:
#' \itemize{
#'  \item \code{std_weights_kDa} Integer with the molecular weight in kDa of the standard proteins.
#'  \item \code{std_elu_fractions} Integer fraction number of the apex of the associated standard protein elution.
#' }
'exampleCalibrationTable'

#' Examplary output of the \code{findComplexFeatures} function
#'
#' Each row in this data.table corresponds to a detected complex co-elution signal. Multiple features per 
#' complex hypotheses can be detected.  
#' @examples 
#'  ## complex feature finding result generation
#'  calibration <- calibrateMW(exampleCalibrationTable)
#'  complexFeatures <- findComplexFeatures(traces = exampleProteinTraces, 
#'                                         complex_hypothesis = exampleComplexHypotheses) 
#'  all.equal(complexFeatures,exampleComplexFeatures)
'exampleComplexFeatures'

#' Examplary output of the \code{findProteinFeatures} function
#'
#' Each row in this data.table corresponds to a detected protein co-elution signal. Multiple features per 
#' protein can be detected.  
#' @examples 
#' ## Protein feature finding result generation
#' calibration <- calibrateMW(exampleCalibrationTable)
#' proteinFeatures <- findProteinFeatures(traces = examplePeptideTracesFiltered) 
#' all.equal(proteinFeatures,exampleProteinFeatures)
'exampleProteinFeatures'

#' Examplary peptide level quantitative PCP-MS dataset in long format.
#' 
#' A data.table containing quantitative peptide level information across multiple
#' fractionations that were analysed in separate MS runs.
#' The structure of this data.table is as follows:
#' \itemize{
#'  \item \code{protein_id} Character string with protein ids (here Uniprot ids).
#'  \item \code{peptide_id} Character string with peptide ids (peptide sequence).
#'  \item \code{filename} Character string with filenames of the separate MS runs.
#'  \item \code{intensity} Intensity values or other quantitative information.
#' }
'examplePCPdataLong'

#' Examplary peptide level quantitative PCP-MS dataset in wide format
#' 
#' A data.table containing quantitative peptide level information across multiple
#' fractionations that were analysed in separate MS runs.
#' The structure of this data.table is as follows:
#' \itemize{
#'  \item \code{protein_id} Character string with protein ids (here Uniprot ids).
#'  \item \code{peptide_id} Character string with peptide ids (peptide sequence).
#'  \item \code{<filename1>} Intensity values or other quantitative information of <filename1>.
#'  \item \code{<filename2>} Intensity values or other quantitative information of <filename2>.
#'  \item \code{...}
#'  \item \code{<filenameX>} Intensity values or other quantitative information of <filenameX>.
#' }
'examplePCPdataWide'

#' Examplary raw peptide level traces object before any annootation
#' 
#' Peptide level traces object as generated directly from the import of a PCP-MS dataset \code{importPCPdata}
#' such as \code{examplePCPdataLong} or \code{examplePCPdataWide}.
#' @examples 
#' ## raw peptide traces object generation from examplary PCP-MS data
#' rawPeptideTraces <- importPCPdata(input_data = examplePCPdataLong,
#'                                   fraction_annotation = exampleFractionAnnotation)
#' all.equal(rawPeptideTraces,examplePeptideTracesUnannotated)
'examplePeptideTracesUnannotated'

#' Examplary peptide level traces object before any filtering
#' #' 
#' Peptide level traces object as generated by annotating a raw peptide level tarces object with 
#' Uniprot annotations and a molecular weight calibration.
#' @examples 
#' ## peptide traces annotation
#' rawPeptideTraces <- examplePeptideTracesUnannotated
#' calibration <- calibrateMW(exampleCalibrationTable)
#' peptideTracesUniprot <- annotateTraces(rawPeptideTraces,
#'                                     trace_annotation = exampleTraceAnnotation)
#' peptideTracesUniprotCalibration <- annotateMolecularWeight(peptideTracesUniprot,
#'                                                                calibration = calibration)
#' all.equal(peptideTracesUniprotCalibration,examplePeptideTraces)
'examplePeptideTraces'

#' Examplary peptide level traces object after filtering
#' 
#' Peptide level traces object as generated by filtering for consecutive id stratches and sibling peptide correlation.
#' @examples 
## filtered peptide traces generation
#' consecutiveFilteredTraces <- filterConsecutiveIdStretches(examplePeptideTraces,
#'                                                           min_stretch_length = 3)
#' sibPepCorrFilteredTraces <- filterBySibPepCorr(consecutiveFilteredTraces,
#'                                                fdr_cutoff = NULL, 
#'                                                absolute_spcCutoff = 0.2)
#' all.equal(sibPepCorrFilteredTraces,examplePeptideTracesFiltered)
'examplePeptideTracesFiltered'

#' Examplary protein level traces object
#' 
#' Protein level traces as, for example, generated by \code{proteinQuantification}.
#' @examples 
#' ## protein quantification
#' proteinTraces <- proteinQuantification(examplePeptideTracesFiltered) 
#' all.equal(proteinTraces,exampleProteinTraces)
'exampleProteinTraces'

#' A sample data.table representing the output format of the OpenSWATH workflow
'exampleOpenSWATHinput'
