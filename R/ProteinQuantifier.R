# For some pretty weird reason
.datatable.aware=TRUE


ProteinQuantifier <- function(pepTraces, No_OF_PEPs = 2 , Keep_Less = FALSE , remove.decoys = TRUE) {
  
  labels <- as.data.frame(pepTraces$traces_annotation)
  labels <- subset(labels, select = c('FullPeptideName','ProteinName'))
  data <- as.data.frame(pepTraces$traces)
  data <- merge(data, labels, by.x = 'id', by.y = 'FullPeptideName',all.x = TRUE)
  data <- as.data.frame(data[,c((2:(ncol(data)-1)),1,ncol(data))])
  
  indx <- sapply(data, is.factor)
  data[indx] <- lapply(data[indx], function(x) as.numeric(as.character(x)))
  
  if (remove.decoys == TRUE) {
    data <- subset(data , substr(data$ProteinName,1,5) != 'DECOY')
  }
  
  proteins <- unique(data$ProteinName)
  nproteins <- length(proteins)
  prot.traces <- data.frame()
  
  for (i in 1:nproteins) {
    
    prot <- proteins[i]
    handle <- subset(data, data$ProteinName == prot)
    handle2 <- subset(handle, select = c('id' , 'ProteinName'))
    handle$ProteinName <- NULL
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
      handle$ProteinName[1] <- prot
      handle <- handle[1,]
      
      prot.traces <- rbind(handle, prot.traces)
      
      rm(handle)
      
    } else if (Keep_Less == TRUE) {
      
      handle$sum <- NULL
      handle <- rbind(colSums(handle[,1:(ncol(handle)-2)]),handle)
      handle$id[1] <- prot
      handle$ProteinName[1] <- prot
      
      handle[1,1:(ncol(handle)-2)] <- (No_OF_PEPs/(ncol(handle)-1)) * handle[1,1:(ncol(handle)-2)]
      
      handle <- handle[1,]
      prot.traces <- rbind(handle, prot.traces)
      
      rm(handle)
    } else {
      rm(handle)
    }
    }
  
  prot.traces$ProteinName <- NULL
  
  traces_annotation <- as.data.frame(pepTraces$traces_annotation)
  traces_annotation <- as.data.table(traces_annotation[(pepTraces$traces_annotation$ProteinName) %in% (prot.traces$id),])
  traces_type <- 'protein'
  fraction_annotation <- as.data.table(pepTraces$fraction_annotation)
  traces <- as.data.table(prot.traces)

  results <- list("traces" = traces,
                  "traces_type" = traces_type,
                  "traces_annotation" = traces_annotation,
                  "fraction_annotation" = fraction_annotation)
  class(results) <- "Traces"
  names(results) <- c("traces", "traces_type", "traces_annotation", "fraction_annotation")
  
  
  return(results)
}

