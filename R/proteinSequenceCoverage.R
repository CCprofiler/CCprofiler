
#' Plot a distribution of the prptides per Protein in your peptideTraces object(s)
#' @param pepTraces_list A list of one or more peptideTraces objects to plot distributions for
#' @param scaled \code{boolean}, Wether to draw a scaled density plot or a histogram
#' @param log \code{boolean}, Wether to plot a logarithmic x-axis
#' @param cutoff \code{numeric}, The limit for the x-axis (peptides per protein)
#' @return A density plot or histogram of the peptides per protein in peptideTraces objects.
#' @details Can be used to compare peptide coverage before and after filtering, or from different searches.
#' Provide a list wit all peptideTraces objects to be compared.
#' @import ggplot2
#' @import data.table
#' @export

plotPeptidesPerProtein <- function(pepTraces_list, scaled = F, log = T, cutoff = 50){
  
  annotation <- lapply(pepTraces_list, "[[","trace_annotation")
  names(annotation) <- names(pepTraces_list)
  pepProtein <- lapply(annotation, function(x) x[,.(peptide_no = .N), by = protein_id])
  pepProtein <- lapply(1:length(pepProtein), function(i) pepProtein[[i]][,name := names(annotation[i])])
  pepProtein <- do.call("rbind", pepProtein)
  if(scaled){
    p <- ggplot(data = pepProtein) +
      geom_density(aes(x = peptide_no, fill = name, colour = name), alpha = 0.3, adjust = 1) 
    
  }else{
    p <- ggplot(data = pepProtein, aes(x = peptide_no, fill = name)) +
      geom_histogram(position = "dodge", binwidth = 1)
  }
  p <- p + scale_x_continuous(limits = c(0, cutoff))
  if(log){
    p <- p + scale_x_log10() 
  }
  
  p
}

#' Calculate the peptide Sequence coverage in every protein by realignment to the protein sequence
#' @param pepTraces peptideTraces object to calculate the sequence coverage from
#' @param proteinFasta Path to a fasta file with Protein Sequences and headers that contain
#' the Protein names in the protein_id column of your peptideTraces object. Usually this will be the
#' fasta file used for library generation in the proteomics search.
#' @return A \code{data.table} containing measured proteins, the length of their AA sequence and the fraction
#' of the sequence covered by peptides
#' @import ggplot2
#' @import data.table
#' @import Biostrings
#' @export


calculateSequenceCoverage <- function(pepTraces, proteinFasta){
  
  protseqs <- readAAStringSet(proteinFasta)
  # proteins <-lapply(pepTraces_list, function(x) unique(x$trace_annotation$protein_id))
  proteins <- unique(pepTraces$trace_annotation$protein_id)
  seq_cov <- sapply(proteins, function(protein){
    #Decoys will not get a match
    ind <- grep(protein, names(protseqs))
    if(length(ind == 1)){
      protseq <- protseqs[[ind]]
      pepseqs <- AAStringSet(pepTraces$trace_annotation[protein_id == protein, gsub("\\(.*?\\)", "", id)])
      #Align the peptides to the protein seq exactly
      match <- matchPDict(pepseqs, protseq)
      # reduce the ranges so overlapping coverage is only counted once
      cov_match <- reduce(unlist(match))
      l <- length(protseq)
      seq_coverage <- sum(width(cov_match)) / l
      
    }else{
      seq_coverage <- NA
      l <- NA
    }
    c(seq_coverage = seq_coverage, length = l)
  })
  seq_cov <- t(seq_cov)
  seq_cov_ <- as.data.table(seq_cov)
  seq_cov_$protein_id <- rownames(seq_cov)
  setcolorder(seq_cov_, c("protein_id","length", "seq_coverage"))
  return(seq_cov_)
}


#' Plot a distribution of Sequence Coverage in your peptideTraces object(s)
#' @param pepTraces_list A list of one or more peptideTraces objects to plot distributions for
#' @param proteinFasta Path to a fasta file with Protein Sequences and headers that contain
#' the Protein names in the protein_id column of your peptideTraces object. Usually this will be the
#' fasta file used for library generation in the proteomics search.
#' @param scaled \code{boolean}, Wether to draw a scaled density plot or a histogram
#' @param log \code{boolean}, Wether to plot a logarithmic x-axis
#' @return A density plot or histogram of the peptides per protein in peptideTraces objects.
#' @details Can be used to compare peptide coverage before and after filtering, or from different searches.
#' Provide a list wit all peptideTraces objects to be compared.
#' @import ggplot2
#' @import data.table
#' @export

plotSequenceCoverage <- function(pepTraces_list, proteinFasta, scaled = TRUE, log = FALSE){
  
  seq_coverage <- lapply(pepTraces_list, calculateSequenceCoverage, proteinFasta)
  seq_coverage <- lapply(1:length(seq_coverage), function(i) seq_coverage[[i]][,Name := names(pepTraces_list)[i]])
  seq_coverage <- do.call("rbind", seq_coverage)
  
  p <- ggplot(data = seq_coverage) 
  if(scale){
    p <- p + geom_density(aes(x = seq_coverage, fill = Name, colour = Name), alpha = 0.3, adjust = 1) 
  }else{
    p <- p + geom_histogram(aes(x = seq_coverage, fill = Name), position = "dodge")
  }
  if(log){
    p <- p + scale_x_log10()
  }
  return(p)
}

#' Estimate the total coverage of all measured proteins, and of all proteins in your fasta database
#' @param sequence_coverage A \code{data.frame} of sequence coverages per protein, as output by 
#' \code{calculateSequenceCoverage}.
#' @param proteinFasta Path to a fasta file with Protein Sequences and headers that contain
#' the Protein names in the protein_id column of your peptideTraces object. Usually this will be the
#' fasta file used for library generation in the proteomics search.
#' @return A list containding the fraction of all AAs covered by peptides of proteins that were measured
#' and of all proteins in the supplied sequence database
#' @import data.table
#' @import Biostrings
#' @export

estimateGlobalSeqenceCoverage <- function(sequence_coverage, proteinFasta){
  protseqs <- readAAStringSet(proteinFasta)
  bases_covered <- sum(sequence_coverage$seq_coverage * sequence_coverage$length, na.rm = T)
  cov_measured <- bases_covered/ sum(sequence_coverage$length, na.rm = T)
  cov_global <- bases_covered / sum(nchar(protseqs))
  return(list(Cov_of_measured_Proteins = cov_measured, Global_Coverage = cov_global))
}


