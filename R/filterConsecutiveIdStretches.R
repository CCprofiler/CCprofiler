# For some pretty weird reason
.datatable.aware=TRUE

#' filterConsecutiveIdStretches filters protein profiling traces based on minimal length of consecutive ID stretches.
#'  
#' @import data.table
#' @param Traces object of class Traces containing the elution profiles that are to be filtered.
#' @param min_stretch_len The minimal length a stretch of contiguous identifications
#'        has to have in order not to be removed.
#' @param remove_empty Logical whether the entries with rowSum == 0 after filtering are removed, defaults to TRUE. 
#' @return An object of class Traces containing the filtered
#'         chromatograms. 
#'
#' @export
#' 

filterConsecutiveIdStretches<-function(Traces, min_stretch_length = 3, remove_empty = TRUE){
  
  #get Traces from container
  
  peptide.traces <- as.data.frame(Traces$peptide.traces)
  
  #Define n as the number of columns
  labels <- peptide.traces[,1:2]
  data <- peptide.traces[,3:ncol(peptide.traces)]
  
  #add 0-column
  data <- cbind(data, dummy = rep(0, nrow(data)))
  
  #Filter
  n <- ncol(data)
  # Set count variable to one
  tmp <- 1
  # going through all rows, for all do:
  for (x in 1:nrow(data)) {
    message(paste("processed", x, "of", nrow(data), "rows"))
    tmp <- 1
    # go through values in all cols, from left to right and ask the following:
    for (i in 1:n){
      #is value in row x, col i = 0?
      if (data[x,i]==0) {
        # if the value is = o, tmp stays = 1
        tmp<-1
      }
      #if value not = 0, then:
      else {
        #look at next value at row x, col i+1, is it 0? If yes, then:
        if (data[x, i+1]==0) {
          # do nothing if the count variable is >3, i.e. there's 4 or more values in a stretch
          if (tmp>(min_stretch_length-1)) {}
          else {
            # replace values if count variable is 3 or smaller at x, i, where the next value is 0. Replace with 0
            for (j in 0:tmp-1) {
              data[x, i-j] <- 0
            }
            tmp <- 1
          }
        }
        # if next value is not = 0, then count tmp+1
        else {
          tmp = tmp+1
        }
      }
    }
  }
  #remove dummy column
  data <- subset(data, select =-dummy)
  #write data to traces object
  peptide.traces.filtered <- cbind(labels, data)
  if (remove_empty) {
    peptide.traces.filtered <- peptide.traces.filtered[rowSums(data) != 0,]
  }
  Traces$peptide.traces <- as.data.table(peptide.traces.filtered)
  return(Traces)
}

# tests
# test<-list(peptide.traces = as.data.table(head(subset(data.wide.pt, select=-decoy), n = 100)))
# 
# filtered <- filterConsecutiveIdStretches(test)
# 
# View(test$peptide.traces)
# View(filtered$peptide.traces)
