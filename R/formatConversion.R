#' Convert a data.table containing traces from wide format to long format. 
#' @param traces.dt A data.table with an id column \code{fragment_id} and
#'        columns of continuously numbered fractions.
#' @return A data.table with columns
#'          \itemize{
#'           \item \code{fragment_id}
#'           \item \code{fraction}
#'           \item \code{intensity}
#'          }
toLongFormat <- function(traces.dt) {
    traces.dt.long <-
        melt(traces.dt, id.var='fragment_id', variable.name='fraction',
             value.name='intensity', variable.factor=FALSE)
    traces.dt.long[, fraction := as.numeric(fraction)]
    traces.dt.long
}


#' Convert a data.table containing traces from long format to wide format. 
#' @param traces.dt A data.table with columns
#'          \itemize{
#'           \item \code{fragment_id}
#'           \item \code{fraction}
#'           \item \code{intensity}
#'          }
#' @return The same data.table casted to wide format.
toWideFormat <- function(traces.dt) {
    dcast(traces.dt, fragment_id ~ fraction, value.var='intensity')
}
