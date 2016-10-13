# For some pretty weird reason
.datatable.aware=TRUE


ProteinQuantifier <- function(pepTraces, No_OF_PEPs = 2 , Keep_Less = FALSE , remove.decoys = TRUE) {
  
  if (remove.decoys) {
    idx_decoys <- grep("^DECOY_",pepTraces$trace_annotation$protein_id)
    pepTraces$traces <- pepTraces$traces[-idx_decoys]
    pepTraces$trace_annotation<- pepTraces$trace_annotation[-idx_decoys]
  }
  
  proteins <- unique(pepTraces$trace_annotation$protein_id)
  nproteins <- length(proteins)
  prot.traces <- data.frame()
  
  data <- as.data.frame(pepTraces$traces)
  data <- merge(data, subset(pepTraces$trace_annotation,select=c("id","protein_id")), by.x = 'id', by.y = 'id',all.x = TRUE)
  data <- as.data.frame(data[,c((2:(ncol(data)-1)),1,ncol(data))])
  
  
  indx <- sapply(data, is.factor)
  data[indx] <- lapply(data[indx], function(x) as.numeric(as.character(x)))
  
  #labels <- as.data.frame(pepTraces$trace_annotation)
  #labels <- subset(labels, select = c('id','protein_id'))
  #data <- as.data.frame(pepTraces$traces)
  #data <- merge(data, labels, by.x = 'id', by.y = 'id',all.x = TRUE)
  #data <- as.data.frame(data[,c((2:(ncol(data)-1)),1,ncol(data))])
  
  #indx <- sapply(data, is.factor)
  #data[indx] <- lapply(data[indx], function(x) as.numeric(as.character(x)))
  
  #if (remove.decoys == TRUE) {
  #  data <- subset(data , substr(data$protein_id,1,5) != 'DECOY')
  #}
  
  #proteins <- unique(data$protein_id)
  #nproteins <- length(proteins)
  #prot.traces <- data.frame()
  
  for (i in 1:nproteins) {
    
    prot <- proteins[i]
    handle <- subset(data, data$protein_id == prot)
    handle2 <- subset(handle, select = c('id' , 'protein_id'))
    handle$protein_id <- NULL
    handle$id <- NULL
    handle <- cbind(handle, handle2)
    rm(handle2)
    
    handle$sum <- as.numeric(rowSums(handle[,1:(ncol(handle)-2)])) 
    handle <- handle[order(-handle$sum),]
    
    if (nrow(handle) >= No_OF_PEPs) {
      
      handle <- handle[1:No_OF_PEPs,]
      handle$sum <- NULL
      
      handle <- rbind(colSums(handle[,1:(ncol(handle)-2)]),handle)
      handle$id[1] <- prot
      handle$protein_id[1] <- prot
      handle <- handle[1,]
      
      prot.traces <- rbind(handle, prot.traces)
      
      rm(handle)
      
    } else if (Keep_Less == TRUE) {
      
      handle$sum <- NULL
      handle <- rbind(colSums(handle[,1:(ncol(handle)-2)]),handle)
      handle$id[1] <- prot
      handle$protein_id[1] <- prot
      
      handle[1,1:(ncol(handle)-2)] <- (No_OF_PEPs/(ncol(handle)-1)) * handle[1,1:(ncol(handle)-2)]
      
      handle <- handle[1,]
      prot.traces <- rbind(handle, prot.traces)
      
      rm(handle)
    } else {
      rm(handle)
    }
    }
  
  prot.traces$protein_id <- NULL
  
  trace_annotation <- data.frame(id=prot.traces$id)
  
  trace_type <- 'protein'
  fraction_annotation <- as.data.table(pepTraces$fraction_annotation)
  traces <- as.data.table(prot.traces)

  results <- list("traces" = traces,
                  "trace_type" = trace_type,
                  "trace_annotation" = trace_annotation,
                  "fraction_annotation" = fraction_annotation)
  class(results) <- "traces"
  names(results) <- c("traces", "trace_type", "trace_annotation", "fraction_annotation")
  
  
  return(results)
}

