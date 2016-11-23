#' Generate binary network from concrete complex hypotheses (e.g. CORUM).
#' @import data.table
#' @param complex_hypothesis data.table in the format of complex hypotheses
#' @return data.table with binary interactions between "a" and "b"
#' @export

generateBinaryNetwork <- function(complex_hypothesis){
  proteins <- unique(complex_hypothesis$protein_id)
  binary_interactions <- data.table(a=character(0),b=character(0))
  for(i in c(1:length(proteins))){
    complexes <- complex_hypothesis$complex_id[which(complex_hypothesis$protein_id == proteins[i])]
    complex_proteins <- unique(complex_hypothesis$protein_id[which(complex_hypothesis$complex_id %in% complexes)])
    complex_proteins <- complex_proteins[complex_proteins != proteins[i]]
    dt <- data.table(a=rep(proteins[i],length(complex_proteins)),b=complex_proteins)
    binary_interactions <- rbind(binary_interactions,dt)
  }
  binary_interactions[]
}

#' Calculate shortest paths between all members of a binary interaction network.
#' This function creates a table that contains all possible combinations of 2 from
#' all input proteins (binomial coeffitient (length(proteins) over 2)).
#' @import data.table
#' @import igraph
#' @param binaryInteractionTable data.table with binary interactions between "a" and "b" as returned from generateBinaryNetwork()
#' @return data.table with three columns (x,y,dist) containing all possible combinations
#' of proteins and their shortest distance as determined by the input binary interaction table
#' @export

calculatePathlength <- function(binaryInteractionTable){
  proteins <- unique(c(binaryInteractionTable$a,binaryInteractionTable$b))
  g <- graph.data.frame(binaryInteractionTable)
  distMatrix <- as.data.frame(shortest.paths(g, v=V(g), to=V(g),weights = NULL, mode="all"))
  combinations <- as.data.table(t(combn(proteins,2)))
  names(combinations) <- c("x","y")
  #combinations <- head(combinations)
  combinations[,rowid := 1:nrow(combinations)]
  combinations[,dist:=distMatrix[x,y],by=rowid]
  combinations[,rowid:=NULL]
  combinations[]
}

#' Generate target hypotheses from distance table. This is done by choosing all
#' max_distance (e.g. 1st) degree interactors for each protein in a distance table.
#' @import data.table
#' @param dist_info data.table with three columns (x,y,dist) containing all possible combinations
#' of proteins and their shortest distance.
#' This is the output from calculatePathlength()
#' @param max_distance numeric maximum distance between to proteins to be considered as being in one complex
#' @return data.table in the format of complex hypotheses
#' @export

generateComplexTargets <- function(dist_info,max_distance=1){
  all_proteins <- unique(c(dist_info$x,dist_info$y))
  complex_table <- data.table(complex_id=character(0),complex_name=character(0),protein_id=character(0))
  unique_complexes <- vector(mode="character")
  t=0
  x=0
  for (i in c(1:length(all_proteins))) {
    print(paste0(i," / ",length(all_proteins)))
    protein <- all_proteins[i]
    combi_sub <- subset(dist_info,(((x == protein) | (y %in% protein)) & (dist<=max_distance)))
    interactors <- sort(unique(c(combi_sub$x,combi_sub$y)))
    # test if exact same complex already exists
    interactor_string <- paste(interactors,collapse ="_")
    if(! interactor_string %in% unique_complexes) {
      x=x+1
      DT <- data.table(complex_id=rep(paste0("c",x),length(interactors)),complex_name=rep(paste0("complex_",x),length(interactors)),protein_id=interactors)
      complex_table <- rbind(complex_table,DT)
      unique_complexes <- c(unique_complexes,interactor_string)
    } else {
      t=t+1
    }
  }
  complex_table[]
}

#' Generate list of proteins for one complex decoy hypothesis
#' @import data.table
#' @param size numeric number of proteins to generate for complex decoy hypothesis
#' @param all_proteins character vector of all possible protein_ids
#' @param dist_info data.table with three columns (x,y,dist) containing all possible combinations
#' of proteins and their shortest distance.
#' This is the output from calculatePathlength()
#' @param min_distance numeric minimum distance between to proteins to be able to be in same complex decoy
#' @return data.table with binary interactions between "a" and "b"
generateDecoys <- function(size,all_proteins,dist_info,min_distance){
  complex_proteins <- vector(mode="character",length=size)
  for(j in c(1:size)){
    combi_sub <- subset(dist_info,(((x %in% complex_proteins) | (y %in% complex_proteins)) & (dist<min_distance)))
    impossible_proteins <- unique(c(combi_sub$x,combi_sub$y))
    possible_proteins <- all_proteins[! all_proteins %in% impossible_proteins]
    complex_proteins[j] <- sample(possible_proteins,1)
  }
  paste(complex_proteins,collapse=";")
}

#' Generate complex decoy hypotheses from target hypotheses and distance table.
#' Selects proteins by min_distance in the network.
#' @import data.table
#' @param target_hypothesis data.table in the format of complex hypotheses
#' @param dist_info data.table with three columns (x,y,dist) containing all possible combinations
#' of proteins and their shortest distance.
#' This is the output from calculatePathlength()
#' @param min_distance numeric minimum distance between to proteins to be able to be in same complex decoy
#' @param seed numeric seed for random number generator
#' @param n.cores numeric number of cores used for parallelization
#' @return data.table with binary interactions between "a" and "b"
#' @export
generateComplexDecoys <- function(target_hypothesis,dist_info,min_distance=2,seed=123,n.cores=1){
  proteins <- unique(target_hypothesis$protein_id)
  target_hypothesis <- target_hypothesis[,list(protein_count=length(complex_id)),by=list(unique.complex_id=complex_id)]
  size=target_hypothesis$protein_count
  cl <- snow::makeCluster(n.cores)
  # setting a seed is absolutely crutial to ensure reproducible results!!!!!!!!!!!!!!!!!!!
  clusterSetRNGStream(cl,seed)
  doSNOW::registerDoSNOW(cl)
  clusterExport(cl, c("generateDecoys"))
  pb <- txtProgressBar(max = length(size), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  #decoy_list <- parLapply(cl=cl,X=size,fun=generateDecoys,all_proteins=proteins,dist_info=dist_info,min_distance=min_distance)
  decoy_list <- foreach(i=seq_along(size),.options.snow = opts) %dopar% {
    x = size[i]
    generateDecoys(size=x,all_proteins=proteins,dist_info=dist_info,min_distance=min_distance)
    }
  close(pb)
  parallel::stopCluster(cl)
  #decoy_list <- lapply(size,generateDecoys,all_proteins=proteins,dist_info=dist_info,min_distance=min_distance)
  decoy_list <- unlist(decoy_list)
  decoy_hypothesis <- data.table(complex_id=paste0("decoy_",seq(1,length(decoy_list),1)),complex_name=paste0("decoy complex ",seq(1,length(decoy_list),1)),protein_id=decoy_list)
  decoy_hypothesis <- decoy_hypothesis[,list(protein = unlist(strsplit(protein_id, ";"))), by=c("complex_id","complex_name")]
  names(decoy_hypothesis)[3] <- "protein_id"
  decoy_hypothesis[]
}