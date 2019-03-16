#' annotateMassDistribution
#' @description Estimates fraction of assembled vs. monomeric mass per id.
#' @param traces An object of type traces. These should be the raw, unprocessed traces.
#' @export
annotateMassDistribution <- function(traces){
  UseMethod("annotateMassDistribution", traces)
}

#' @describeIn annotateMassDistribution Estimates fraction of assembled vs. monomeric mass per id.
#' @export
annotateMassDistribution.traces <- function(traces){
  .tracesTest(traces)

  if(! ("protein_mw" %in% names(traces$trace_annotation))) {
    stop("annotateMassDistribution requires information on the monomer molecular weight. Please see annotateTraces.")
  }

  if(! ("molecular_weight" %in% names(traces$fraction_annotation))) {
    stop("annotateMassDistribution requires information on the molecular weight assigned to each fraction. Please see annotateMolecularWeight.")
  }

  intMat <- getIntensityMatrix(traces)
  mw_info <- traces$trace_annotation
  max_trace <- max(traces$fraction_annotation$id)
  mw_info[,assembly_boundary := 2*protein_mw]
  mw_info[is.na(assembly_boundary), assembly_boundary := 0]
  mw_info[,assembly_boundary_fraction := unlist(lapply(assembly_boundary,function(x){max(traces$fraction_annotation[molecular_weight>=x]$id)}))]

  # Get intensity sum
  mw_info[,sum_assembled := sum(intMat[id,1:assembly_boundary_fraction-1]), id]
  mw_info[,sum_monomeric := sum(intMat[id,assembly_boundary_fraction:max_trace]), id]
  mw_info[,sum := sum(intMat[id,]), id]
  mw_info[is.na(assembly_boundary), sum_assembled := 0]
  mw_info[is.na(assembly_boundary), sum_monomeric := 0]
  mw_info[, sum_assembled_norm := 1/sum*sum_assembled]
  mw_info[, sum_monomeric_norm := 1/sum*sum_monomeric]
  #mw_info[, test := sum_monomeric_norm+sum_assembled_norm]

  traces$trace_annotation <- mw_info
  .tracesTest(traces)
  return(traces)
}

#' @describeIn annotateMassDistribution Estimates fraction of assembled vs. monomeric mass per id.
#' @export
annotateMassDistribution.tracesList <- function(tracesList){
  .tracesListTest(tracesList)
  tracesListRes <- lapply(tracesList, annotateMassDistribution.traces)
  class(tracesListRes) <- "tracesList"
  .tracesListTest(tracesListRes)
  return(tracesListRes)
}

#' getMassAssemblyChange
#' @description Estimates change of assembled vs. monomeric mass per id.
#' @param tracesList An object of type traces.list including assembled mass
#' annotation as produced by \code{annotateMassDistribution}.
#' @param design_matrix data.table, design matrix describing the architecture of the tracesList object.
#' @param plot logical if plot should be created, default is FALSE.
#' @param name character string specifying the name of output if PDF=TRUE, default is "beta_pvalue_histogram".
#' @param PDF logical if PDF should be created, default is FALSE.
#' @export
getMassAssemblyChange <- function(tracesList, design_matrix,
                                  compare_between = "Condition",
                                  quantLevel = "protein_id",
                                  plot = FALSE,
                                  PDF = FALSE,
                                  name = "beta_pvalue_histogram"){
  .tracesListTest(tracesList)
  samples <- unique(design_matrix$Sample)
  if(! all(samples %in% names(tracesList))) {
    stop("tracesList and design_matrix do not match. Pleas check sample names.")
  }
  if (! "sum_assembled_norm" %in% names(tracesList[[1]]$trace_annotation)) {
    stop("No assembled mass annotation available, please run annotateMassDistribution first.")
  }
  if (! quantLevel %in% names(tracesList[[1]]$trace_annotation)) {
    stop("quantLevel not available in provided traces.")
  }

  res <- lapply(names(tracesList), function(tr){
    vals <- subset(tracesList[[tr]]$trace_annotation,select=c(quantLevel,"sum_assembled_norm"))
    vals[,Sample := tr]
    return(vals)
  })

  res <- do.call(rbind, res)
  
  if(length(unique(design_matrix$Replicate)) < 2) {
    if (quantLevel == "protein_id") {
      res_cast <- dcast(res, formula = protein_id ~ Sample, value.var=c("sum_assembled_norm"))
    } else if (quantLevel == "proteoform_id") {
      res_cast <- dcast(res, formula = proteoform_id ~ Sample, value.var=c("sum_assembled_norm"))
    } else {
      stop("Functionality only available for quantLevel proetin_id or proteoform_id.")
    }
    #res_cast[, change := log2(get(samples[1])/(get(samples[2])))]
    #res_cast[change=="NaN", change := 0]
    res_cast[, medianDiff := get(samples[1])-(get(samples[2]))]
    res_cast[, betaPval := 1]
    res_cast[, betaPval_BHadj := 1]
    res_cast[, testOrder := paste0(samples[1],".vs.",samples[2])]
    res_cast <- subset(res_cast, select = c("protein_id","medianDiff","betaPval", "betaPval_BHadj","testOrder"))
    return(res_cast[])
  } else {
    res <- merge(res, design_matrix, by.x="Sample", by.y="Sample_name")
    res[,n_conditions:=length(unique(Condition)), by=c("protein_id")]
    res[,n:=.N, by=c("protein_id")]
    res[,replicates_perCondition:=.N, by=c("protein_id", "Condition")]
    #res[,sum_assembled_norm := ifelse(sum_assembled_norm>0.999,sum_assembled_norm-0.001,sum_assembled_norm)]
    #res[,sum_assembled_norm := ifelse(sum_assembled_norm<0.001,sum_assembled_norm+0.001,sum_assembled_norm)]
    res[,sum_assembled_norm_t := (sum_assembled_norm * (n - 1) + 0.5)/n, by=c("protein_id")]
    res[,unique_perCondition := length(unique(round(sum_assembled_norm, digits = 3))), by=c("protein_id","Condition")]
    
    diff <- res[, { 
      samples = unique(.SD[,get(compare_between)])
      median = .SD[, .(m = median(sum_assembled_norm)), by = .(get(compare_between))]
      medianDiff = median[get==samples[1]]$m - median[get==samples[2]]$m
      n_conditions <- unique(.SD$n_conditions)
      n_perCondition <- min(.SD$replicates_perCondition)
      n_unique_perCondition <- min(.SD$unique_perCondition)
      if( (n_conditions > 1) & (n_perCondition > 1) & (n_unique_perCondition > 1) ) {
        model = betareg(.SD$sum_assembled_norm_t ~ .SD$Condition)
        stat = lrtest(model)
        p = stat$`Pr(>Chisq)`[2]
      } else {
        p = 2
      }
      .(medianDiff = medianDiff, betaPval = p, testOrder = paste0(samples[1],".vs.",samples[2]))}, 
      by = .(get(quantLevel))]
    
    diff[betaPval==2, betaPval := NA ]
    diff[, betaPval_BHadj := p.adjust(betaPval, method = "fdr")]
    
    if(plot==TRUE){
      if(PDF){
        pdf(paste0(name,".pdf"))
      }
      hist(diff$betaPval, breaks = 100)
      if(PDF){
        dev.off()
      }
    }
    
    setnames(diff, "get(quantLevel)", quantLevel)
    tests <- subset(diff, select = c("protein_id","medianDiff","betaPval", "betaPval_BHadj","testOrder"))
    
    return(tests[])
    
  }
}

#' plotMassAssemblyChange between 2 conditions
#' @param assamblyTest data.table, a data.table with test statistics.
#' An assamblyTest can be produced with \code{getMassAssemblyChange}.
#' @param change_cutoff Numeric fold change cutoff, default is 0.2.
#' @param name character string specifying the name of output if PDF=TRUE, default is "massAssemblyChange".
#' @param PDF logical if PDF should be created, default is FALSE.
#' @return plot
#' @export
plotMassAssemblyChange <- function(assamblyTest, change_cutoff=0.2, name="massAssemblyChange", PDF=FALSE){
  if (PDF) {
    pdf(paste0(name,".pdf"))
  }
  no_change <- nrow(subset(assamblyTest, abs(medianDiff) <= change_cutoff))
  more_assembled <- nrow(subset(assamblyTest, medianDiff < -change_cutoff))
  less_assembled <- nrow(subset(assamblyTest, medianDiff > change_cutoff))
  pie(c(no_change,more_assembled,less_assembled),
      labels=c(paste("no assembly change \n abs(change) < ",change_cutoff,"\n",no_change),
               paste0("more assembled \n change < -",change_cutoff,"\n",more_assembled),
               paste0("less assembled \n change > ",change_cutoff,"\n",less_assembled)),
      main = paste0(unique(assamblyTest$testOrder)))

  h <- ggplot(assamblyTest, aes(x=medianDiff)) + geom_histogram(binwidth = 0.05) + theme_classic()
  print(h)
  if (PDF) {
    dev.off()
  }
}
