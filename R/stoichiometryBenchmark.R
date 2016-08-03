stoichiometryBenchmark.plot <- function(swf){
  up <- fread("/Volumes/ibludau-1/SEC/stoichiometry/uniprot-human_032016.csv")
  
  corum_st <- fread("/Volumes/ibludau-1/SEC/stoichiometry/corum_complexes_with_stoichs.tsv")
  proteingroups <- strsplit(corum_st$protein_stoichiometry, split = ",")
  
  masses <- list()
  totalmasses <- c()
  for (i in 1:length(proteingroups)){
    proteins <- proteingroups[[i]]
    merged <- merge(data.table(Entry = proteins), up, by = "Entry", all.x = TRUE)
    masses[[i]] <- gsub(",","", merged$Mass)
    masses[[i]] <- as.numeric(masses[[i]])/1000
    totalmasses[i] <- sum(masses[[i]])
  }
  corum_st[, masses:=sapply(masses, function(x){paste(x, collapse = ",")})]
  corum_st[, cumulative_mass_PDBstoichiometries:=totalmasses]
  
  ### results from peakpicking
  best.features <- getBestComplexFeature(swf)
  detected.features <- subsetComplexFeatures(best.features,min_completeness = 1)
  detected.features[,max.stoichiometry_estimated := NA]
  for(i in 1:nrow(detected.features)){
    max <- max(as.integer(unlist(strsplit(detected.features$stoichiometry_estimated[[i]], split = ";"))))
    detected.features$max.stoichiometry_estimated[i] <- max
  }
  detected.features <- subset(detected.features,max.stoichiometry_estimated <= 2)
  
  # Merge with the complex feature table and plot MS-estimated-MWs vs PDB-calculated-MWs
  features_and_mws_completeobs <- merge(corum_st,detected.features, by = "complex_id")
  
  features_and_mws_completeobs <- unique(features_and_mws_completeobs,by="complex_name.x")
  
  mod1 = lm(cumulative_mass_PDBstoichiometries~complex_mw_estimated, data = features_and_mws_completeobs)
  modsum = summary(mod1)
  r2 = modsum$adj.r.squared
  mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
  
  q <- ggplot(features_and_mws_completeobs, aes(x = cumulative_mass_PDBstoichiometries,
                                                y = complex_mw_estimated,colour=as.factor(n_subunits_detected)))
  q + geom_point() +
    scale_colour_brewer(palette = "Set1",name="number of subunits") +
    xlab("MW from PDB-extracted stoichiometries [kDa]") +
    ylab("MW from MS estimated stoichiometries [kDa]") +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    ggtitle("Stoichiometry Estimates") +
    ylim(0,500) +
    xlim(0,500) +
    geom_text(aes(label=complex_name.x), size=2,hjust=-0.1)
}

