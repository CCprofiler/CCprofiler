#' createGroupsFromNetworkTable create protein groups around network nodes 
#' with all proteins connected by 1 edge (first degree interactors)
#' @import data.table
#' @import doSNOW
#' @param interaction_table Interaction network in edge table format (data.frame)
#' @param id_column_1 Name of the column containing protein identifier #1
#' @param id_column_2 Name of the column containing protein identifier #2
#' @return A data.frame containing group name (seed node name) and its first line interactors
#' e.g. for use as "complex hypotheses" for  
#' @export

createGroupsFromNetworkTable <- function(interaction_table, n.cores = detectCores()-1 , id_column_1 = "UniprotA", id_column_2 = "UniprotB"){
  # helper function to filter for the first line interactors of one protein
  
  interaction_table <- as.data.table(interaction_table)
  
  filterInteractors<- function(protein_id, network_table, prot_a_colid = id_column_1, prot_b_colid = id_column_2){
    
	# uniprot<-gene_uniprot_map[gene_uniprot_map[,1] %in% genename,2]
	a_indices <- which(t(subset(network_table, select = prot_a_colid)) %in% protein_id)
	b_indices <- which(t(subset(network_table, select = prot_b_colid)) %in% protein_id)
	
	prim_interactions <- network_table[c(a_indices, b_indices),]
	
	colnum_a <- which(names(network_table) == prot_a_colid)
	colnum_b <- which(names(network_table) == prot_b_colid)
	
	prim_interactors <- unique(c(as.character(t(prim_interactions[,colnum_a, with = FALSE])),
                                 as.character(t(prim_interactions[,colnum_b, with = FALSE]))))
	
	return(prim_interactors)
  }
  
  ## Apply to all protein_ids in a parallelized fashion
  # 1. set up compute nodes
  cl <- snow::makeCluster(n.cores)
  doSNOW::registerDoSNOW(cl)
  
  # Apply filterInteractors using foreach
  nodes <- unique(c(t(interaction_table[,id_column_1,with=FALSE])), t(interaction_table[,id_column_2,with=FALSE])) # nodes = IA network nodes
  n.nodes <- length(nodes)
  
  modules.n1 <- foreach (i = seq_along(nodes), .combine=c, .export = "data.table") %dopar% {
    interactors <- filterInteractors(nodes[i], interaction_table)
    list(data.table(complex_id = rep(i, length(interactors)),
                    complex_name = rep(paste0(nodes[i], "_interactors_degree1"), length(interactors)),
                    protein_id = interactors))
  }
  
  # collect results in data.frame
  return(rbindlist(modules.n1))
 
  }

  