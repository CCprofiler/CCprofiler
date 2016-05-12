# For some pretty weird reason
.datatable.aware=TRUE

#' filterConsecutiveIdStretches filters protein profiling traces based on
#' minimal length of consecutive ID stretches.
#'  
#' @import data.table
#' @param traces object of class traces containing the elution profiles that
#'     are to be filtered.
#' @param min_stretch_len The minimal length a stretch of contiguous
#'     identifications has to have in order not to be removed.
#' @param remove_empty Logical whether the entries with rowSum == 0 after
#'     filtering are removed, defaults to TRUE. 
#' @return An object of class traces containing the filtered
#'     chromatograms. 
#' @export
filterConsecutiveIdStretches<-function(traces,
                                       min_stretch_length=3,
                                       remove_empty=TRUE) {
  # Get traces from container
  peptide.traces <- as.data.frame(traces$peptide.traces)
  # Define n as the number of columns
  labels <- peptide.traces[, 1:2]
  data <- peptide.traces[, 3:ncol(peptide.traces)]
  # Add 0-column
  data <- cbind(data, dummy=rep(0, nrow(data)))
  # Filter
  n <- ncol(data)
  # Set count variable to one
  tmp <- 1
  # Going through all rows, for all do:
  for (x in 1:nrow(data)) {
    message(paste('processed', x, 'of', nrow(data), 'rows'))
    tmp <- 1
    # Go through values in all cols, from left to right and ask the following:
    for (i in 1:n) {
      # Is value in row x, col i = 0?
      if (data[x,i] == 0) {
        # If the value is = o, tmp stays = 1
        tmp<-1
      }
      # If value not = 0, then:
      else {
        # Look at next value at row x, col i+1, is it 0? If yes, then:
        if (data[x, i+1] == 0) {
          # Do nothing if the count variable is >3, i.e. there's 4 or more values in a stretch
          if (tmp > (min_stretch_length-1)) {}
          else {
            # Replace values if count variable is 3 or smaller at x, i, where the next value is 0. Replace with 0
            for (j in 0:tmp-1) {
              data[x, i-j] <- 0
            }
            tmp <- 1
          }
        }
        # If next value is not = 0, then count tmp+1
        else {
          tmp <- tmp+1
        }
      }
    }
  }

  # Remove dummy column
  data <- subset(data, select =-dummy)
  # Write data to traces object
  peptide.traces.filtered <- cbind(labels, data)
  if (remove_empty) {
    peptide.traces.filtered <- peptide.traces.filtered[rowSums(data) != 0, ]
  }
  traces$peptide.traces <- as.data.table(peptide.traces.filtered)
  return(traces)
}

# tests
# test<-list(peptide.traces = as.data.table(head(subset(data.wide.pt, select=-decoy), n = 100)))
# 
# filtered <- filterConsecutiveIdStretches(test)
# 
# View(test$peptide.traces)
# View(filtered$peptide.traces)
