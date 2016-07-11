#' Get a matrix of fragment intensity values. 
#' @param traces.obj An object of type \code{FragmentTraces}.
#' @return A matrix of intensity values.
#' @export
getIntensityMatrix <- function(traces.obj) {
    ids <- traces.obj$traces$id
    intensity.mat <- as.matrix(subset(traces.obj$traces,
                                      select=-id))
    rownames(intensity.mat) <- ids
    intensity.mat
}
