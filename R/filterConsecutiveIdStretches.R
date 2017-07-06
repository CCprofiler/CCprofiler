# Due to: http://stackoverflow.com/questions/24501245/data-table-throws-object-not-found-error
.datatable.aware=TRUE

#' Filter Consecutive Id Stretches in Chromaograms
#' @description Filter PCP traces based on a
#' minimal length of consecutive ID stretches.
#' @import data.table
#' @param traces An object of class traces.
#' @param min_stretch_length Numeric integer, the minimal length a stretch of continuous
#'     identifications has to have in order not to be removed.
#' @param remove_empty Logical, whether the entries with rowSum == 0 after
#'     filtering are removed, defaults to TRUE.
#' @return An object of class traces containing the filtered
#'     chromatograms.
#' @export
#' @example
#' 
## Load example data
# tracesRaw <- examplePeptideTraces
# 
# tracesFilteredConSec <- filterConsecutiveIdStretches(tracesRaw,
#                                                      min_stretch_length = 3,
#                                                      remove_empty = T)
filterConsecutiveIdStretches<-function(traces,
                                       min_stretch_length=3,
                                       remove_empty=TRUE) {
  ## Test traces
  .tracesTest(traces)
  ## Get traces from container
  peptideTraces <- getIntensityMatrix(traces)
  
  id <- rownames(peptideTraces)
  intensity <- peptideTraces
  ## Add 0-column
  intensity <- cbind(intensity, dummy=rep(0, nrow(intensity)))
  
  ## Define n as the number of columns
  n <- ncol(intensity)
  ## Set count variable to one
  tmp <- 1
  ## Going through all rows, for all do:
  for (x in 1:nrow(intensity)) {
    #message(paste('processed', x, 'of', nrow(intensity), 'peptides'))
    tmp <- 1
    ## Go through values in all cols, from left to right and ask the following:
    for (i in 1:n) {
      ## Is value in row x, col i = 0?
      if (intensity[x,i] == 0) {
        ## If the value is = 0, tmp is set to 1
        tmp <- 1
      }else {# If value not = 0, then:
        ## Look at next value at row x, col i+1, is it 0? 
        if (intensity[x, i+1] == 0) { #If yes, then:
          ## Do nothing if the count variable is > min_stretch_length
          if (tmp <= (min_stretch_length-1)){ # if count variable is min_stretch_length or smaller at x
            ## Replace values where the next value is 0. Replace with 0
            for (j in 0:tmp-1) {
              intensity[x, i-j] <- 0
            }
            tmp <- 1
          }
        }else {# If next value is not = 0, then count tmp+1
          tmp <- tmp + 1
        }
      }
    }
  }

  # Remove dummy column
  intensity <- subset(intensity, select =-dummy)

  # Write intensity to traces object
  peptideTracesFiltered <- as.data.table(intensity)
  peptideTracesFiltered[,id := rownames(intensity)]

  tracesAnnotation <- traces$trace_annotation

  if (remove_empty) {
    idx <- which(rowSums(intensity) != 0)
    peptideTracesFiltered <- peptideTracesFiltered[idx]
    tracesAnnotation <- tracesAnnotation[idx]
  }

  traces$traces <- peptideTracesFiltered
  traces$trace_annotation <- tracesAnnotation
  .tracesTest(traces)
  return(traces)
}
