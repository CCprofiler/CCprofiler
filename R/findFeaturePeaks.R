
#' A helper function to extend a list of complex features with additional
#' information.
#'
#' @param features A data.table of complex feature candidates with the
#'        following format:
#'        \itemize{
#'         \item \code{subgroup} A semicolon-separated list of protein
#'                               identifiers.
#'         \item \code{left_sec} The left boundary of the feature.
#'         \item \code{right_sec} The right boundary of the feature.
#'         \item \code{score} The intra-feature correlation.
#'        }
#' @param trace.mat A matrix where rows correspond to protein traces.
#'        This is the matrix that was used to find complex features.
#' @param protein.names A character vector specifying the protein identifiers
#'        belonging to the rows in \code{trace.mat}.
#' @param protein.mw.conc A data.table that stores the molecular weight and
#'        estimate of the absolute abundance for each subunit.
#'        \itemize{
#'         \item \code{protein_id}
#'         \item \code{protein_mw}
#'         \item \code{protein_concentration}
#'        }
#' @param protein.abundances A numeric vector holding estimates of the absolute
#'        abundance of each subunit.
#' @return The same data.table as the input argument extended with the
#'         following columns:
#'         \itemize{
#'          \item \code{n_subunits} The number of subunits in the feature.
#'          \item \code{stoichiometry} The intensity-based stoichiometry.
#'          \item \code{mw_estimated} The estimated molecular weight.
#'          \item \code{mw_apparent} The apparent mw.
#'          \item \code{mw_delta} The delta mw.
#'         }
findFeaturePeaks <- function(features, trace.mat,protein.names,protein.mw.conc) {
    # Compute the number of subunits in each complex feature
    features[, n_subunits := length(strsplit(as.character(subgroup), ';')[[1]]),by=subgroup]

    traces.dt <- data.table(trace.mat)
    setnames(traces.dt,as.character(seq(1,ncol(traces.dt),1)))
    traces.dt[, protein_id := protein.names]

    # create complex trace by summing up intensities
    library('pracma')  # @TODO put into dependencies
    features.new <- lapply(seq(1:nrow(features)), function(i){
      #print(i)
       feature=features[i]
       subunits <- strsplit(feature$subgroup, ';')[[1]]
       # Extract only the traces for the subunits that make up this feature.
       subunit.traces <- subset(traces.dt, protein_id %in% subunits)
       subunit.traces[,protein_id:=NULL]
       # create complex.trace by summing traces for all subunits
       complex.trace <- subunit.traces[, lapply(.SD, sum, na.rm=FALSE)]
       complex.trace.mat <- as.matrix(complex.trace)
       #smoothing
       complex.trace.SG <-savgol(complex.trace.mat, fl=11, forder = 4, dorder = 0)
      # peak picking
       #complex.peaks <- data.table(findpeaks(complex.trace.mat[1,],minpeakdistance=5,nups=3,ndowns=3))
       complex.peaks <- findpeaks(complex.trace.SG,minpeakdistance=5,nups=3,ndowns=3)
       if (is.null(dim(complex.peaks))){
         complex.peaks <- data.table(intensity=complex.peaks[1],
           apex=complex.peaks[2],
           left=complex.peaks[3],
           right=complex.peaks[4])
       } else {
         complex.peaks <- data.table(complex.peaks)
         setnames(complex.peaks,c("intensity","apex","left","right"))
       }

       # select peaks within boundaries of correlation based window
       sel_peaks <- which((complex.peaks$left>=feature$left_sec) & (complex.peaks$right<=feature$right_sec))
       if (length(sel_peaks > 0)) {
         complex.peaks <- complex.peaks[sel_peaks]
         # only peak with highest intensity
         complex.peak <- complex.peaks[which(complex.peaks$intensity==max(complex.peaks$intensity)),]
         complex.peak$area <- sum(complex.trace.mat[1,complex.peak$left:complex.peak$right])
       } else {
         apex=feature$left_sec+((feature$right_sec - feature$left_sec)/2)
         apex=unique(c(floor(apex),ceiling(apex)))
         intensity=mean(complex.trace.mat[apex])
         complex.peak <- data.table(intensity=intensity,apex=mean(apex),left=feature$left_sec,right=feature$right_sec)
         complex.peak$area <- sum(complex.trace.mat[1,complex.peak$left:complex.peak$right])
       }
       complex.peak
     })

     features[, intensity := unlist(lapply(features.new, function(x) x$intensity))]
     features[, apex := unlist(lapply(features.new, function(x) x$apex))]
     features[, left := unlist(lapply(features.new, function(x) x$left))]
     features[, right := unlist(lapply(features.new, function(x) x$right))]
     features[, area := unlist(lapply(features.new, function(x) x$area))]

     # Estimate the molecular weight of the complex by using information
     # gathered by calibrating the column with molecules of known molecular
     # weight. This is known as the 'apparent' molecular weight.
     features[, mw_apparent := convertSECToMW(left+((right - left)/2)),
              by=subgroup] # @CORRECTED for central sec fraction

     # Produce a long list version of the trace matrix since its more convenient
     # for downstream processing.
     traces.long <- melt(traces.dt, id.vars='protein_id', value.name='intensity',
                         variable.name='fraction', variable.factor=FALSE)
     traces.long[, fraction := as.numeric(fraction)]
     setkey(traces.long, protein_id)

     # Next, create a data.table that holds information on subunits.
     protein.info <- traces.long[, list(total_intensity=sum(intensity)), by=protein_id]
     # Merge the data.table holding the molecular weight and abs. abundance
     # estimates and the data.table with the total intensity of each protein.
     setkey(protein.mw.conc, protein_id)
     setkey(protein.info, protein_id)
     protein.info <- protein.info[protein.mw.conc]
     setkey(protein.info, protein_id)

     # Estimate the stoichiometry of each complex feature
     # as well as its molecular mass.
     res <- lapply(seq(1:nrow(features)), function(i){
        feature=features[i]
        subunits <- strsplit(feature$subgroup, ';')[[1]]
        # Extract only the traces for the subunits that make up this feature.
        subunit.traces <- traces.long[subunits]
        # Extract the relevant information from the protein.info DT.
        subunit.mws <- protein.info[subunits, protein_mw]
        subunit.abundances <- protein.info[subunits, protein_concentration]
        subunit.total.intensities <- protein.info[subunits, total_intensity]
        # Call the helper function to estimate the molecular mass of the
        # complex and its stoichiometry.
        res <- estimateComplexMass(subunit.traces,
          feature$left,
          feature$right,
          subunit.total.intensities,
          subunit.mws,
          subunit.abundances)
     })

     features[, mw_estimated := unlist(lapply(res, function(x) x$mw_estimated))]
     features[, stoichiometry := unlist(lapply(res, function(x) x$stoichiometry))]

     # Produce a score that measures how well both estimates agree.
     features[, mw_delta := abs(mw_estimated - mw_apparent)]

     print(features)  ##### @ERROR THIS IS A BUG IN R I GUESS!!!!!
     return(features)
   }
