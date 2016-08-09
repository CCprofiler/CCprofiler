#setwd("/Volumes/ibludau-1/SEC/SECprofiler")
#require(devtools)
#load_all()
#document()
#protein_ids = c("P20248","P24941","P18858","P09874","P09884","P28340","Q07864","P35251","P35250","P27694","P15927","P35244","P11387")

SECexplorer_processing <- function(protein_ids){
    #@TODO ID mapping  
    protein.traces <- traces.obj
    protein.traces <- subset(protein.traces,trace_ids = protein_ids)
    complex.table <- data.table(complex_id=rep(1,length(protein_ids)),complex_name=rep("1",length(protein_ids)),protein_id=protein_ids)
    swf <- findComplexFeatures(traces.obj = protein.traces,
                               complex.protein.assoc = complex.table,
                               corr.cutoff = 0.9,
                               window.size = 15,
                               parallelized = FALSE)
    
    res = resultsToTable(swf)
    list(traces=protein.traces,table=res)
}
