#' Get a matrix of fragment intensity values. 
#' @param trace.obj An object of type \code{FragmentTraces}.
#' @return A matrix of intensity values.
#' @export
getIntensityMatrix <- function(trace.obj) {
    ids <- trace.obj$traces$id
    intensity.mat <- as.matrix(subset(trace.obj$traces,
                                      select=-id))
    rownames(intensity.mat) <- ids
    intensity.mat
}
