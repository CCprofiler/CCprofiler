#' generateBinaryNetwork
#' @description Generate binary network from concrete complex hypotheses (e.g. CORUM).
#' @import data.table
#' @param complex_hypothesis data.table in the format of complex hypotheses
#' @return data.table with binary interactions between "a" and "b"
#' @export
generateBinaryNetwork <- function(complex_hypothesis){
  complexes <- unique(complex_hypothesis$complex_id)
  binary_interactions <- data.table(a=character(0),b=character(0))
  # get all binary interactions for each complex in the complex hypothesis
  for(i in c(1:length(complexes))){
    binary_interactions <- appendBinaryInteractions(complex=complexes[i],hypothesis=complex_hypothesis,binary_interactions=binary_interactions)
  }
  setkey(binary_interactions,NULL)
  unique(binary_interactions) # remove duplicates
  binary_interactions <- unique(binary_interactions)
  binary_interactions[]
}

#' appendBinaryInteractions
#' @description Append binary interaction table with additional complexes
#' @import data.table
#' @param complex character string with complex_id
#' @param hypothesis data.table in the format of complex hypotheses
#' @param data.table with binary interactions between "a" and "b"
#' @return data.table with binary interactions between "a" and "b" appended by current complex interactions
appendBinaryInteractions <- function(complex,hypothesis,binary_interactions){
  complex_proteins <- unique(hypothesis$protein_id[which(hypothesis$complex_id == complex)]) # get all proteins within complex
  binaries <- as.list(data.table(combn(complex_proteins,2))) # create all binary combinations of the proteins in complex
  binaries_sort <- lapply(binaries,FUN=sort) # sort alphanumerically to make it possible to remove duplicates at the end
  binaries_sort_table <- data.table(t(as.data.table(binaries_sort)))
  colnames(binaries_sort_table) <- c("a","b")
  binary_interactions <- rbind(binary_interactions,binaries_sort_table)
  binary_interactions
}

#' calculatePathlength
#' @description Calculate shortest paths between all members of a binary interaction network.
#' This function creates a table that contains all possible combinations of 2 from
#' all input proteins (binomial coeffitient (length(proteins) over 2)).
#' @import data.table
#' @import igraph
#' @param binaryInteractionTable data.table with binary interactions between "a" and "b" as returned from generateBinaryNetwork()
#' @return data.table with three columns (x,y,dist) containing all possible combinations
#' of proteins and their shortest distance as determined by the input binary interaction table
#' @export
calculatePathlength <- function(binaryInteractionTable){
  proteins <- unique(c(binaryInteractionTable$a,binaryInteractionTable$b)) # get all protein ids in interaction table
  g <- graph_from_data_frame(binaryInteractionTable, directed = FALSE) # create undirected graph using the igraph package
  distMatrix <- distances(g) #as.data.frame(shortest.paths(g, v=V(g), to=V(g),weights = NULL, mode="all")) # calculate shortest path between all vertices
  combinations <- melt(distMatrix)
  combinations <- as.data.table(combinations)
  names(combinations) = c("x","y","dist")
  combinations$x <- as.character(combinations$x)
  combinations$y <- as.character(combinations$y)
  combinations[]
}

#' generateComplexTargets
#' @description Generate target hypotheses from distance table. This is done by choosing all
#' max_distance (e.g. 1st) degree interactors for each protein in a distance table.
#' @import data.table
#' @param dist_info data.table with three columns (x,y,dist) containing all possible combinations
#' of proteins and their shortest distance.
#' This is the output from calculatePathlength()
#' @param max_distance numeric maximum distance between to proteins to be considered as being in one complex
#' @param redundancy_cutoff numeric maximum overlap distance between two hypotheses (0=identical,1=subset,between 1 and 2=some shared subunits, 2=no shared subunits), default=1
#' @return data.table in the format of complex hypotheses
#' @export
generateComplexTargets <- function(dist_info,max_distance=1,redundancy_cutoff=1){
  all_proteins <- unique(c(dist_info$x,dist_info$y))
  #complex_table <- data.table(complex_id=character(0),complex_name=character(0),protein_id=character(0))
  initial_complexes <- data.table(complex_id=character(0),subunits_detected=character(0),n_subunits=character(0))
  for (i in c(1:length(all_proteins))) {
    #print(paste0(i," / ",length(all_proteins)))
    protein <- all_proteins[i]
    combi_sub <- subset(dist_info,(((x == protein) | (y == protein)) & (dist<=max_distance)))
    interactors <- sort(unique(c(combi_sub$x,combi_sub$y)))
    interactor_string <- paste(interactors,collapse =";")
    DT <- data.table(complex_id=protein,subunits_detected=interactor_string,n_subunits=length(interactors))
    initial_complexes <- rbind(initial_complexes,DT)
  }
  # remove redundant hypotheses
  dist_hyp <- getDistanceMatrix(initial_complexes)
  clust_hyp <- complexClustering(initial_complexes,dist_hyp)
  tree_cut=cutree(clust_hyp,h=redundancy_cutoff)
  initial_complexes[,consecutive_feature_identifier := .I]
  initial_complexes[,unique_feature_identifier := 0]
  unique_feature_identifier=0
  unique_tree_groups <- unique(tree_cut)
  for (tree_group in unique_tree_groups) {
    unique_feature_identifier = unique_feature_identifier+1
    tree_complex_ids <- which(tree_cut==tree_group)
    initial_complexes$unique_feature_identifier[tree_complex_ids] = unique_feature_identifier
  }
  initial_complexes <- initial_complexes[order(unique_feature_identifier,-n_subunits)]
  initial_complexes[,complex_name := paste(complex_id,collapse="_"),by="unique_feature_identifier"]
  initial_complexes_unique <- unique(initial_complexes,by="unique_feature_identifier")
  initial_complexes_unique[,complex_id := paste0("c",.I)]
  initial_complexes_unique[,complex_name := paste0("complex_",complex_name)]
  complex_table <- subset(initial_complexes_unique,select=c("complex_id","complex_name","subunits_detected"))
  complex_table <- complex_table[,list(protein_id = unlist(strsplit(subunits_detected, ";"))), by=c("complex_id","complex_name")]
  complex_table[]
}

#' collapseHypothesis
#' @description Remove redundancy in existing complex hypotheses
#' @import data.table
#' @param hypothesis data.table with complex hypotheses
#' @param redundancy_cutoff numeric maximum overlap distance between two hypotheses (0=identical,1=subset,between 1 and 2=some shared subunits, 2=no shared subunits), default=1
#' @return data.table in the format of complex hypotheses
#' @export
collapseHypothesis <- function(hypothesis,redundancy_cutoff=1){
  hypothesis[,n_subunits := length(protein_id),by="complex_id"]
  hypothesis[,subunits_detected := paste(protein_id,collapse=";"),by="complex_id"]
  hypothesis <- unique(hypothesis,by="complex_id")
  hypothesis[,protein_id := NULL]
  dist_hyp <- getDistanceMatrix(hypothesis)
  clust_hyp <- complexClustering(hypothesis,dist_hyp)
  tree_cut=cutree(clust_hyp,h=redundancy_cutoff)
  hypothesis[,consecutive_feature_identifier := .I]
  hypothesis[,unique_feature_identifier := 0]
  unique_feature_identifier=0
  tree_features <- names(tree_cut)
  unique_tree_groups <- unique(tree_cut)
  for (tree_group in unique_tree_groups) {
    unique_feature_identifier = unique_feature_identifier+1
    tree_complex_ids <- which(tree_cut==tree_group)
    hypothesis$unique_feature_identifier[tree_complex_ids] = unique_feature_identifier
  }
  hypothesis <- hypothesis[order(unique_feature_identifier,-n_subunits)]
  hypothesis[,complex_name := paste(complex_name,collapse="; "),by="unique_feature_identifier"]
  hypothesis_unique <- unique(hypothesis,by="unique_feature_identifier")
  hypothesis_unique <- subset(hypothesis_unique,select=c("complex_id","complex_name","subunits_detected"))
  hypothesis_unique <- hypothesis_unique[,list(protein_id = unlist(strsplit(subunits_detected, ";"))), by=c("complex_id","complex_name")]
  hypothesis_unique
}



#' generateDecoys
#' @description Generate list of proteins for one complex decoy hypothesis
#' @import data.table
#' @param size numeric number of proteins to generate for complex decoy hypothesis
#' @param all_proteins character vector of all possible protein_ids
#' @param dist_info data.table with three columns (x,y,dist) containing all possible combinations
#' of proteins and their shortest distance.
#' This is the output from calculatePathlength()
#' @param min_distance numeric minimum distance between to proteins to be able to be in same complex decoy
#' @return data.table with binary interactions between "a" and "b"
#' @export
generateDecoys <- function(size,all_proteins,dist_info,min_distance){
  complex_proteins <- vector(mode="character",length=size)
  for(j in c(1:size)){
    combi_sub <- subset(dist_info,(((x %in% complex_proteins) | (y %in% complex_proteins)) & (dist<min_distance)))
    impossible_proteins <- unique(c(combi_sub$x,combi_sub$y))
    possible_proteins <- all_proteins[! all_proteins %in% impossible_proteins]
    if (length(possible_proteins) > 0){
      complex_proteins[j] <- sample(possible_proteins,1)
    } else {
      stop(paste0("Network connectivity is too high for decoy generation at a min_distance of ",min_distance, ". Please select a lower min_distance cutoff."))
    }
  }
  paste(complex_proteins,collapse=";")
}

#' generateComplexDecoys
#' @description Generate complex decoy hypotheses from target hypotheses and distance table.
#' Selects proteins by min_distance in the network.
#' @import data.table
#' @param target_hypothesis data.table in the format of complex hypotheses
#' @param dist_info data.table with three columns (x,y,dist) containing all possible combinations
#' of proteins and their shortest distance.
#' This is the output from calculatePathlength()
#' @param min_distance numeric minimum distance between to proteins to be able to be in same complex decoy
#' @param append logical if to append decoys to the input target hypotheses, default is FALSE
#' @param seed numeric seed for random number generator
#' @param parallelized logical, if the computation should be done in parallel, default=FALSE
#' @param n_cores numeric number of cores used for parallelization
#' @return data.table with binary interactions between "a" and "b"
#' @export
generateComplexDecoys <- function(target_hypothesis,dist_info,min_distance=2,append=FALSE,seed=123,parallelized=FALSE,n_cores=1){
  input = target_hypothesis
  proteins <- unique(target_hypothesis$protein_id)
  target_hypothesis <- target_hypothesis[,list(protein_count=length(complex_id)),by=list(unique.complex_id=complex_id)]
  size=target_hypothesis$protein_count
  if (parallelized) {
    cl <- snow::makeCluster(n_cores)
    # setting a seed is absolutely crutial to ensure reproducible results!!!!!!!!!!!!!!!!!!!
    clusterSetRNGStream(cl,seed)
    doSNOW::registerDoSNOW(cl)
    clusterExport(cl, c("generateDecoys"))
    pb <- txtProgressBar(max = length(size), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    set.seed(seed)
    #decoy_list <- parLapply(cl=cl,X=size,fun=generateDecoys,all_proteins=proteins,dist_info=dist_info,min_distance=min_distance)
    decoy_list <- foreach(i=seq_along(size),.options.snow = opts) %dopar% {
      x = size[i]
      generateDecoys(size=x,all_proteins=proteins,dist_info=dist_info,min_distance=min_distance)
      }
    close(pb)
    parallel::stopCluster(cl)
  } else {
    set.seed(seed)
    pb <- txtProgressBar(max = length(size), style = 3)
    decoy_list <- foreach(i=seq_along(size)) %do% {
      setTxtProgressBar(pb, i)
      x = size[i]
      generateDecoys(size=x,all_proteins=proteins,dist_info=dist_info,min_distance=min_distance)
      }
    close(pb)
  }
  #decoy_list <- lapply(size,generateDecoys,all_proteins=proteins,dist_info=dist_info,min_distance=min_distance)
  decoy_list <- unlist(decoy_list)
  decoy_hypothesis <- data.table(complex_id=paste0("DECOY_",seq(1,length(decoy_list),1)),complex_name=paste0("DECOY_",seq(1,length(decoy_list),1)),protein_id=decoy_list)
  decoy_hypothesis <- decoy_hypothesis[,list(protein = unlist(strsplit(protein_id, ";"))), by=c("complex_id","complex_name")]
  names(decoy_hypothesis)[3] <- "protein_id"
  if (append == TRUE) {
    appended_hypotheses <- rbind(input,decoy_hypothesis)
    appended_hypotheses[]
  } else if (append == FALSE) {
    decoy_hypothesis[]
  }  else {
    stop("Invalid parameter for append. Please select TRUE or FALSE.")
  }
}
