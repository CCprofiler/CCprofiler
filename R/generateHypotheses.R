#' 
#' Generate a binary network from complexes.
#' @description Generate binary network from concrete complex hypotheses (e.g. CORUM).
#' @import data.table
#' @param complex_hypotheses data.table specifying complex hypotheses.
#' Must have the following columns:
#' \itemize{
#' \item complex_id: character strings, a unique id for every complex
#' \item protein_id: character strings, the protein id, e.g. Uniprot id
#' }
#' @return data.table with binary interactions between "a" and "b"
#' @export
#' @examples 
#' ## Load example Data
#' complexHypotheses <- exampleComplexHypotheses
#' 
#' ## Generate the binary network
#' binaryInteractions <- generateBinaryNetwork(complex_hypotheses = complexHypotheses)
#' ##Inspect the result
#' binaryInteractions
#' 
#' 
generateBinaryNetwork <- function(complex_hypotheses){
  setkey(complex_hypotheses, "complex_id")
  complexes <- unique(complex_hypotheses$complex_id)
  ## Get every combination of binary interactions within a complex
  binary_interactions <- do.call("rbind", lapply(complexes, function(x){
    t(combn(complex_hypotheses[complex_id == x, protein_id], m = 2))
  }))
  ## Keep only unique interactions
  binary_interactions <- t(apply(binary_interactions,1, sort))
  binary_interactions <- as.data.table(binary_interactions)
  binary_interactions <- unique(binary_interactions)
  
  names(binary_interactions) <- c("a", "b")
  setkey(binary_interactions,NULL)
  
  return(binary_interactions)
}


#' Calculate pathlength from binary interactions
#' @description Calculate shortest paths between all members of a binary interaction network.
#' This function creates a table that contains all possible combinations of 2 from
#' all input proteins (binomial coefficient (length(proteins) over 2)).
#' @import data.table
#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph distances
#' @param binaryInteractionTable data.table with binary interactions between "a" and "b" as returned from generateBinaryNetwork()
#' @return data.table with three columns containing all possible combinations of proteins:
#' \itemize{
#'   \item x,y: Protein ids of the input proteins (all possible combinations)
#'   \item dist: distance (shortest path between interactors)
#'         as determined by the input binary interaction table
#' }
#' @export
#' @examples 
#' ## Load example data
#' complexHypotheses <- exampleComplexHypotheses
#' 
#' ## Generate the binary network
#' binaryInteractions <- generateBinaryNetwork(complex_hypotheses = complexHypotheses)
#' 
#' ## Calculate the path lengths
#' shortestPaths <- calculatePathlength(binaryInteractions)
#' 

calculatePathlength <- function(binaryInteractionTable){
  ## get all protein ids in interaction table
  proteins <- unique(c(binaryInteractionTable$a,binaryInteractionTable$b)) 
  ## create undirected graph using the igraph package
  g <- graph_from_data_frame(binaryInteractionTable, directed = FALSE)
  ## calculate shortest path between all vertices
  distMatrix <- distances(g) 
  ## Output as a long format data.table
  combinations <- melt(distMatrix)
  combinations <- as.data.table(combinations)
  names(combinations) = c("x","y","dist")
  combinations$x <- as.character(combinations$x)
  combinations$y <- as.character(combinations$y)
  
  return(combinations)
}

#' Generate complex target-hypotheses
#' @description Generate target hypotheses from distance table. This is done by choosing all
#' max_distance (e.g. 1st) degree interactors for each protein in a distance table.
#' @import data.table
#' @param dist_info data.table with three columns containing all possible combinations of proteins:
#' \itemize{
#'   \item x,y: Protein ids of the input proteins (all possible combinations)
#'   \item dist: distance as determined by the input binary interaction table
#' }
#' This is the output from \code{\link[SECprofiler]{calculatePathlength}}
#' @param max_distance Numeric, maximum distance between to proteins to be considered
#' as being in one complex. Defaults to 1.
#' @param redundancy_cutoff Numeric, maximum overlap distance between two hypotheses
#' (0=identical,1=subset,between 1 and 2=some shared subunits, 2=no shared subunits).
#' Defaults to 1.
#' 
#' @return data.table in the format of complex hypotheses.
#' Has the following columns:
#' \itemize{
#' \item complex_id: character strings, a unique id for every complex
#' \item protein_id: character strings, the protein id, e.g. Uniprot id
#' }
#' 
#' @export
#' @examples
#' ## Load example Data
#' complexHypotheses <- exampleComplexHypotheses
#' 
#' ## Generate the binary network
#' binaryInteractions <- generateBinaryNetwork(complex_hypotheses = complexHypotheses)
#' 
#' ## Calculate the path lengths
#' shortestPaths <- calculatePathlength(binaryInteractions)
#' 
#' ## Generate a set of target hypotheses from the binary interactions
#' targetHypotheses_ <- generateComplexTargets(dist_info = shortestPaths,
#'                                             max_distance = 1,
#'                                             redundancy_cutoff = 1)
#' 

generateComplexTargets <- function(dist_info,
                                   max_distance=1,
                                   redundancy_cutoff=1){
  
  all_proteins <- unique(c(dist_info$x,dist_info$y))
  
  ## Extract one complex hypotheis for every protein contained
  initial_complexes <- lapply(all_proteins, function(protein){
    
    combi_sub <- subset(dist_info,(((x == protein) | (y == protein)) & (dist<=max_distance)))
    interactors <- sort(unique(c(combi_sub$x,combi_sub$y)))
    interactorString <- paste(interactors,collapse =";")
    data.table(complex_id=protein,
               subunits_detected=interactorString,
               n_subunits=length(interactors))
  })
  initial_complexes <- do.call("rbind", initial_complexes)
  
  ## Remove redundant hypotheses
  complex_table <- .collapseWideHypothesis(hypothesis = initial_complexes,
                                           redundancy_cutoff = redundancy_cutoff)
  return(complex_table)
}

#' Collapse redundant hypotheses
#' @description Remove redundancy in existing complex hypotheses
#' @import data.table
#' @param hypothesis data.table with complex hypotheses
#' Must have the following columns:
#' \itemize{
#' \item complex_id: character strings, a unique id for every complex
#' \item protein_id: character strings, the protein id, e.g. Uniprot id
#' }
#' Additionaly a complex_name column is allowed.
#' @param redundancy_cutoff Numeric, maximum overlap distance between two hypotheses
#' \itemize{
#'   \item 0 = identical 
#'   \item 1 = subset
#'   \item between 1 and 2 = some shared subunits
#'   \item 2 = no shared subunits 
#' }
#' Defaults to 1.
#' @return data.table in the format of complex hypotheses
#' Has the following columns:
#' \itemize{
#'   \item complex_id: character strings, a unique id for every complex
#'   \item complex_name: character strings, Name of the complex
#'   \item protein_id: character strings, the protein id, e.g. Uniprot id
#' }
#' @export
#' @examples 
#' ## load example data
#' complexHypotheses <- exampleComplexHypotheses
#' ## Collapse redundancies
#' complexHypothesesCollapsed <- collapseHypothesis(complexHypotheses)
#' 
#' 
collapseHypothesis <- function(hypothesis,
                               redundancy_cutoff=1){
  hyp <- copy(hypothesis)
  hyp[,n_subunits := length(protein_id), by="complex_id"]
  hyp[,subunits_detected := paste(protein_id,collapse=";"),by="complex_id"]
  hyp <- unique(hyp, by="complex_id")
  hyp[,protein_id := NULL]
  
  hypothesis_unique <- .collapseWideHypothesis(hyp, redundancy_cutoff = redundancy_cutoff)
  return(hypothesis_unique)
}



#' Generate one comlpex decoy-hypothesis
#' @description Generate list of proteins for one complex decoy hypothesis
#' @import data.table
#' @param size Numeric, number of proteins to generate for complex decoy hypothesis
#' @param all_proteins Character, vector of all possible protein_ids
#' @param dist_info data.table with three columns (x,y,dist) containing all possible combinations
#' of proteins and their shortest distance.
#' This is the output from calculatePathlength()
#' @param min_distance Numeric, minimum distance between to proteins to be
#' able to be in same complex decoy. Defaults to 1
#' @param n_tries Numeric integer, Number of attempts to generate a decoy hypothesis until the
#' algorithm quits and reports an error. Higher values allow for decoy generation with high-connectivity
#' target networks, but take longer to compute.
#' @return character string with a complex hypothesis, protein ids are separated by ';'
#' 
generateDecoys <- function(size,
                           all_proteins,
                           dist_info,
                           min_distance = 1,
                           n_tries = 3){
  success <- F
  n <- 1
  while(success == F & n <= n_tries){
    complex_proteins <- vector(mode="character",length=size)
    
    for(j in c(1:size)){
      combi_sub <- dist_info[(x %in% complex_proteins) | (y %in% complex_proteins)]
      impossible_proteins <- unique(c(combi_sub$x,combi_sub$y))
      possible_proteins <- all_proteins[! all_proteins %in% impossible_proteins]
      if (length(possible_proteins) > 0){
        complex_proteins[j] <- sample(possible_proteins,1)
        success <- T
      }else{
        success <- F
        #   stop(paste0("Network connectivity is too high for decoy generation at a min_distance of ",
        #               min_distance, ". Please select a lower min_distance cutoff."))
      }
    }
    n <- n + 1
  }
  if(success){
    return(paste(complex_proteins,collapse=";"))
  }else{
    stop(paste0("Network connectivity is too high for decoy generation at a min_distance of ",
                min_distance, ". Please select a lower min_distance cutoff."))
  }
  
}

#' Generate comlpex decoy-hypotheses
#' @description Generate complex decoy hypotheses from target hypotheses and distance table.
#' Selects proteins by min_distance in the network.
#' @import data.table
#' @param target_hypotheses data.table in the format of complex hypotheses
#' Must have the following columns:
#' \itemize{
#' \item complex_id: character strings, a unique id for every complex
#' \item protein_id: character strings, the protein id, e.g. Uniprot id
#' }
#' @param dist_info data.table with three columns containing all possible combinations of proteins:
#' \itemize{
#'   \item x,y: Protein ids of the input proteins (all possible combinations)
#'   \item dist: distance as determined by the input binary interaction table
#' }
#' This is the output from \code{\link[SECprofiler]{calculatePathlength}}
#' @param min_distance Numeric, minimum distance between to proteins to be
#' able to be in same complex decoy.
#' @param append Logical, wether to append decoys to the input target hypotheses, default is FALSE.
#' @param seed Numeric, seed for random number generator.
#' @param parallelized Logical, if the computation should be done in parallel, default=FALSE.
#' @param n_cores Numeric, number of cores used for parallelization.
#' @param n_tries Numeric integer, Number of attempts to generate a decoy hypothesis until the
#' algorithm quits and reports an error. Higher values allow for decoy generation with high-connectivity
#' target networks, but take longer to compute.
#' @return data.table in the form of complex hypotheses.
#' Has the following columns:
#' \itemize{
#' \item complex_id: character strings, a unique id for every complex
#' \item complex_name: name of the complex (in this case identical with complex id)
#' \item protein_id: character string, the protein id, e.g. Uniprot id
#' }
#' @export
#' @examples 
#' ## Load example Data
#' complexHypotheses <- exampleComplexHypotheses
#' ## Generate the binary network
#' binaryInteractions <- generateBinaryNetwork(complex_hypotheses = complexHypotheses)
#' ## Calculate the path lengths
#' shortestPaths <- calculatePathlength(binaryInteractions)
#' 
#' ## Generate the Decoys
#' decoys <- generateComplexDecoys(target_hypotheses = complexHypotheses,
#'                                 dist_info = shortestPaths,
#'                                 min_distance = 2)
#' ## Inspect the resulting decoys
#' decoys
#' 

generateComplexDecoys <- function(target_hypotheses,
                                  dist_info,
                                  min_distance=2,
                                  append=FALSE,
                                  seed=123,
                                  parallelized=FALSE,
                                  n_cores=1,
                                  n_tries = 3){
  
  input <- copy(target_hypotheses)
  proteins <- unique(target_hypotheses$protein_id)
  target_hypotheses <- target_hypotheses[,.(protein_count=length(complex_id)),
                                         by=.(unique.complex_id=complex_id)]
  size <- target_hypotheses$protein_count
  ## Reduce the dist_info size to increase speed
  dist_info_subs <- subset(dist_info, dist < min_distance)
  setkeyv(dist_info_subs, c("x", "y"))
  if (parallelized) {
    cl <- snow::makeCluster(n_cores)
    # setting a seed is absolutely crutial to ensure reproducible results!
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
      generateDecoys(size = x, all_proteins = proteins,
                     dist_info = dist_info_subs,
                     min_distance = min_distance,
                     n_tries = n_tries)
    }
    close(pb)
    parallel::stopCluster(cl)
  } else {
    set.seed(seed)
    pb <- txtProgressBar(max = length(size), style = 3)
    decoy_list <- foreach(i=seq_along(size)) %do% {
      setTxtProgressBar(pb, i)
      x = size[i]
      generateDecoys(size = x, all_proteins = proteins,
                     dist_info = dist_info_subs,
                     min_distance = min_distance,
                     n_tries = n_tries)
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
