#' #' #' Cyclic loess normalization
#' @description Normalize peptide traces across replicates, conditions and fractions
#' @import data.table
#' @import preprocessCore
#' @import limma
#' @import evobiR
#' @import plyr
#' @param traces_list traces list object of type peptide
#' @param window integer size of slding window, default = 3
#' @param step integer steps for sliding window, default = 1
#' @param plot logical, default = TRUE
#' @param PDF logical, defalt = FALSE
#' @param name character strimg specifying pdf name, default = "normalizeByCyclicLoess"
#' @export
#'

normalizeByCyclicLoess <- function(traces_list, window = 3, step = 1, plot = TRUE, PDF = TRUE, name = "normalizeByCyclicLoess") {
  .tracesListTest(traces_list, type = "peptide")
  trace_intensities_long <- lapply(traces_list, extractvaluesForNorm)
  combi_table <- rbindlist(trace_intensities_long, use.names=TRUE, fill=FALSE, idcol="sample")
  combi_table[, filename := paste0(sample,"_",fraction_number)]
  combi_table$fraction_number <- as.numeric(combi_table$fraction_number)
  combi_table <- unique(combi_table)

  combi_table_forPlot <- copy(combi_table)
  combi_table_forPlot[, total_intensity:=sum(intensity), by=c("filename","sample")]
  combi_table_forPlot <- unique(subset(combi_table_forPlot, select=c("fraction_number","total_intensity","sample")))
  pnormdata<-ggplot(combi_table_forPlot, aes(x=fraction_number, y=total_intensity, group=sample)) +
    geom_line(aes(color=sample)) +
    geom_point(aes(color=sample)) +
    theme_classic()
  ggsave(pnormdata,filename=paste0(name,"_priorNormalization.pdf"),width=7,height=3.5)

  combi_table[, intensity := log2(intensity)]
  combi_table[, intensity := ifelse(intensity == -Inf, NA, intensity)]
  #combi_table[, intensity := ifelse(intensity < 0.000001, NA, intensity)]
  combi_table_norm <- normalize_sn(combi_table, window, step)
  combi_table_norm[, intensity := 2^(intensity)]
  combi_table_norm[, intensity := ifelse(intensity == 1, 0, intensity)]

  saveRDS(combi_table_norm,"combi_table_norm.rds")

  combi_table_toMerge <- subset(combi_table, select = c("filename","id", "sample", "fraction_number"))
  combi_table_norm_final <- merge(combi_table_toMerge, combi_table_norm, by=c("filename","id"))

  combi_table_norm_forPlot <- copy(combi_table_norm_final)
  combi_table_norm_forPlot[, total_intensity:=sum(intensity), by=c("filename","sample")]
  combi_table_norm_forPlot <- unique(subset(combi_table_norm_forPlot, select=c("fraction_number","total_intensity","sample")))
  pnormdata<-ggplot(combi_table_norm_forPlot, aes(x=fraction_number, y=total_intensity, group=sample)) +
    geom_line(aes(color=sample)) +
    geom_point(aes(color=sample)) +
    theme_classic()
  ggsave(pnormdata,filename=paste0(name,"_postNormalization.pdf"),width=7,height=3.5)

  list_norm <- split(combi_table_norm_final, by="sample")
  list_norm_wide <- lapply(list_norm, dcast_backToTraces)

  traces_list_norm <- copy(traces_list)
  sample_names <- names(traces_list)
  for(s_name in sample_names){
    traces_list_norm[[s_name]]$traces <- list_norm_wide[[s_name]]
    traces_list_norm[[s_name]]$fraction_annotation <- subset(traces_list_norm[[s_name]]$fraction_annotation, id %in% names(traces_list_norm[[s_name]]$traces))
    traces_list_norm[[s_name]]$trace_annotation <- subset(traces_list_norm[[s_name]]$trace_annotation, id %in% traces_list_norm[[s_name]]$traces$id)
  }
  .tracesListTest(traces_list_norm, type = "peptide")
  return(traces_list_norm)
}

dcast_backToTraces <- function(normData){
  normData_sub <- unique(subset(normData, select = c("id","fraction_number","intensity")))
  normData_wide <- data.table::dcast(normData_sub, id ~ fraction_number, value.var = "intensity", drop = TRUE)
  if (ncol(normData_wide) < (max(unique(normData_sub$fraction_number)))+1) {
    missing <- seq(1,max(unique(normData_sub$fraction_number)),1)[which(! seq(1,max(unique(normData_sub$fraction_number)),1) %in% names(normData_wide))]
    for (m in missing){
      normData_wide[, as.character(eval(m)) := NA]
    }
  }
  for (j in seq_len(ncol(normData_wide))){
    set(normData_wide,which(is.na(normData_wide[[j]])),j,0)
  }
  setcolorder(normData_wide, c(seq(1,(ncol(normData_wide)-1),1),"id"))
  setkey(normData_wide, "id")
  normData_wide$id <- as.character(normData_wide$id)
  return(normData_wide)
}

extractvaluesForNorm <- function(traces){
  intensities_long <- data.table::melt(traces$traces, var.id = "id", variable.name = "fraction_number", value.name = "intensity")
  return(intensities_long)
}

normalize_sn <- function(X, window, step) {
  mx<-dcast(X, id~filename, value.var='intensity', sum)
  mx[mx<0.00001]=NA

  mxs<-as.matrix(mx[,-1])
  rownames(mxs)<-mx$id

  id_mapping<-unique(X[,c("filename","fraction_number")])

  max_sec <- max(X$fraction_number)
  windows_sets<-SlidingWindow("data.frame",c(0:max_sec+1), window, step)

  #lmxn<-lapply(windows_sets,function(X){normalizeMedianValues(mxs[,subset(id_mapping, fraction_number %in% X)$filename])})
  lmxn<-lapply(windows_sets,function(X){normalizeCyclicLoess(mxs[,subset(id_mapping, fraction_number %in% X)$filename])})

  lln<-do.call("rbind",lapply(lmxn, melt, na.rm=TRUE))
  names(lln)<-c("id", "filename", "intensity")

  lln_dt <- as.data.table(lln)
  lln_dt[,mean_intensity := mean(intensity, na.rm=T), by=c("id","filename")]
  lln_dt_sub <- unique(subset(lln_dt, select = c("id","filename","mean_intensity")))
  names(lln_dt_sub)<-c("id", "filename", "intensity")
  #lxn<-ddply(lln, .(id,filename),function(X){mean(X$intensity)})
  #names(lxn)<-c("id", "filename", "intensity")
  #return(lxn)
  return(lln_dt_sub)
}
