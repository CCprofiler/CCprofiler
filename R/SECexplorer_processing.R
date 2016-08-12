#setwd("/Volumes/ibludau-1/SEC/SECprofiler")
#require(devtools)
#load_all()
#document()
#protein_ids = c("P20248","P24941","P18858","P09874","P09884","P28340","Q07864","P35251","P35250","P27694","P15927","P35244","P11387")
#protein_ids = c("P61201","Q9UNS2","Q9BT78","Q92905","Q7L5N1","Q99627","Q13098","Q9UBW8","P61201","Q9UNS2","Q9BT78","Q92905","Q7L5N1","Q99627","Q13098","Q9H9Q2")

SECexplorer_processing <- function(protein_ids,traces.obj){
    #@TODO ID mapping?
    traces.obj <- subset(traces.obj,trace_ids = protein_ids)
    complex.table <- data.table(complex_id=rep(1,length(protein_ids)),complex_name=rep("1",length(protein_ids)),protein_id=protein_ids)
    swf <- findComplexFeatures(traces.obj = traces.obj,
                               complex.protein.assoc = complex.table,
                               corr.cutoff = 0.95,
                               window.size = 6,
                               parallelized = FALSE)
    res = resultsToTable(swf)
    list(traces=traces.obj,features=res)
}

#res = SECexplorer_processing(protein_ids = protein_ids)

#for(i in 1:nrow(res$features)){
#  plot.complexFeatures(res$features[i],res$traces,complexID = 1,monomerMW=FALSE)
#}
