#' Get a matrix of fragment intensity values. 
#' @param traces.obj An object of type \code{FragmentTraces}.
#' @return A matrix of intensity values.
#' @export
getIntensityMatrix <- function(traces.obj) {
    ids <- traces.obj$traces$id
    intensity.mat <- as.matrix(sapply(subset(traces.obj$traces,
                                      select=-id),as.numeric))
    rownames(intensity.mat) <- ids
    intensity.mat
}
