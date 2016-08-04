#' Plot the result of the sliding window algorithm.
#' @param sw.result An object of type 'complexFeaturesSW'.
#' @param trace.mat The same matrix that was passend to
#'        \code{findComplexFeatures}.
#' @param protein.names The names of proteins whose traces are contained in
#'        the matrix \code{trace.mat}.
#' @param n.largest Only plot the n largest subgroups. Optional.
#' @export
plot.complexFeatures <- function(res,
                                 traces.obj,
                                 complexID,
                                 monomerMW=FALSE,
                                 log=FALSE) {


    features <- subset(res, complex_id == complexID)
    proteins <- unique(unlist(strsplit(features$subunits_annotated, split = ";")))
    traces <- subset(traces.obj,trace_ids = proteins)
    traces.long <- toLongFormat(traces$traces)
    
    detectedSubunits = strsplit(features$subunits_detected, ';')[[1]]
    sel_inComplex <- which(traces.long$id %in% detectedSubunits)
    traces.long$inComplex <- "FALSE"
    traces.long$inComplex[sel_inComplex] <- "TRUE"
    traces.long$inComplex = factor(traces.long$inComplex,levels=c("TRUE","FALSE"))
    # add monomer MW
    subunitMW = as.numeric(strsplit(features$monomer_sec, ';')[[1]])
    subunitMW.dt = data.table(id=detectedSubunits,mw=subunitMW)
    
    traces.long = merge(traces.long,subunitMW.dt,by="id",all.x = TRUE)
    
    setkey(traces.long, inComplex, id)
    
    complexName = unique(features$complex_name)[1]
    n_annotatedSubunits = features$n_subunits_annotated[1]
    n_signalSubunits = features$n_subunits_with_signal[1]
    n_detectedSubunits = features$n_subunits_detected[1]
    complexCompleteness = round(features$completeness[1],digits=2)
    
    
    p <- ggplot(traces.long) +
      geom_line(aes_string(x='fraction', y='intensity', color='id',linetype='inComplex')) +
      ggtitle(paste0(complexName,
                     "\nAnnotated subunits: ",n_annotatedSubunits,
                     "   Subunits with signal: ",n_signalSubunits,
                     "\nCoeluting subunits: ",n_detectedSubunits,
                     "   Completeness: ",complexCompleteness)) +
      xlab('fraction') +
      ylab('intensity')
    
    if (log) {
        p <- p + scale_y_log10('log(intensity)')
    }

#    p + geom_vline(aes(xintercept=left_sec, color=subgroup),
#                   data=found.features) +
#        geom_vline(aes(xintercept=right_sec, color=subgroup),
#                   data=found.features) +
#        theme(legend.position='none')


    p <- p + geom_rect(data=features,aes(xmin = left_pp, xmax = right_pp, ymin = -Inf, ymax = Inf, fill=subunits_detected),alpha = 0.25)
    p <- p + geom_vline(data=features,aes(xintercept=left_sw), colour="grey",linetype="longdash")
    p <- p + geom_vline(data=features,aes(xintercept=right_sw), colour="grey",linetype="longdash")
    p <- p + geom_vline(data=features,aes(xintercept=apex), colour="black",linetype="solid")
    p <- p + geom_vline(data=features,aes(xintercept=complex_sec_estimated), colour="black",linetype="longdash")
    if (monomerMW==TRUE){
      p <- p + geom_vline(data=subset(traces.long,!is.na(mw)),aes(xintercept=mw,colour=id),linetype="solid")
    }
    print(p)
}
