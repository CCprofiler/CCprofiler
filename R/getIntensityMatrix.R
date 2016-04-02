#' Get a matrix of fragment intensity values. 
#' @param trace.obj An object of type \code{FragmentTraces}.
#' @return A matrix of intensity values.
#' @export
getIntensityMatrix <- function(trace.obj) {
    fragment.ids <- trace.obj$traces$fragment_id
    intensity.mat <- as.matrix(subset(trace.obj$traces,
                                      select=-fragment_id))
    rownames(intensity.mat) <- fragment.ids
    intensity.mat
}
