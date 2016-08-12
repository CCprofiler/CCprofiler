#' Some protein traces of e4 that were filtered. TESTING DATA.
#'protein.traces'

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
'corum.complex.protein.assoc'

#' An example object of class \code{FragmentTraces}.
'traces.obj'

#' A data.table that holds stores for each protein its molecular weight as
#' annotated in UniProt as well as an estimate on its absolute abundance within
#' the sample.
#' Its structure is as follows:
#' \itemize{
#'  \item \code{protein_id} The Uniprot identifier of the protein.
#'  \item \code{protein_mw} The molecular weight.
#'  \item \code{protein_concentration} The absolute concentration.
#' }
'protein.mw.conc'
