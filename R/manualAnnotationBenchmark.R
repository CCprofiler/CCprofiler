manualAnnotationBenchmark.plot <- function(swf){
  #load("/Volumes/ibludau-1/SEC/manual.annotations.final.rda")
  load("/Volumes/ibludau-1/SEC/corum_db/manual.annotations.new.rda")
  manual.annotations.final <- manual.annotations.new
  
  true.complexes <- unique(manual.annotations.final$complex_id)
  
  best.features <- getBestComplexFeature(swf)
  detected.features <- subsetComplexFeatures(best.features,min_completeness = 0.5)
  detected.complexes <- detected.features$complex_id
  
  #input.complexes = unique(corum.complex.protein.assoc$complex_id)
  input.complexes = unique(filtered_corum_table$complex_id)
  
  rates <- estimate_TPR_and_FPR(detected.complexes = detected.complexes,true.complexes = true.complexes,input.complexes = input.complexes)
  TPR=rates$TPR
  FPR=rates$FPR
  
  rates_df <- data.table(completeness = seq(0,1,0.05), TPR=NA, FPR=NA)
  for(i in 1:nrow(rates_df)){
    detected.features <- subsetComplexFeatures(best.features,min_completeness = rates_df$completeness[i])
    detected.complexes <- detected.features$complex_id
    rates <- estimate_TPR_and_FPR(detected.complexes = detected.complexes,true.complexes = true.complexes,input.complexes = input.complexes)
    rates_df$TPR[i] <- rates$TPR
    rates_df$FPR[i] <- rates$FPR
  }
  #plot(x=rates_df$FPR,y=rates_df$TPR)
  
  pl <- ggplot(data=rates_df,aes(x=FPR,y=TPR,colour=completeness)) +
    geom_point() +
    ylim(0,1) +
    xlim(0,1) +
    scale_colour_gradient(low="green", high="blue")
  print(pl)
  
  
  #res = data.table(TPR=TPR, FPR=FPR)
  #write.table(res,"/IMSB/ra/ibludau/SEC/benchmark_15_90_newQuant_LooseBoundaries_noPeak_2nd.txt",sep="\t",quote=F,row.names=F,col.names=T)
  
  TP_complexes <- detected.complexes[which(detected.complexes %in% true.complexes)]
  RT_df <- data.frame(detected = vector(mode="numeric", length=length(TP_complexes)),true = vector(mode="numeric", length=length(TP_complexes)))
  
  i=0
  for(complex in TP_complexes) {
    i=i+1
    detected_features <- subset(detected.features,complex_id == complex)
    RT_detected = detected_features$apex
    
    true_features <- subset(manual.annotations.final,complex_id == complex)
    sel_closest_RT <- which(abs(true_features$rt-RT_detected) == min(abs(true_features$rt-RT_detected)))[1]
    RT_true <- true_features$rt[sel_closest_RT]
    
    RT_df$detected[i] = RT_detected
    RT_df$true[i] = RT_true
  }
  
  mod1 = lm(true~detected, data = RT_df)
  modsum = summary(mod1)
  r2 = modsum$adj.r.squared
  mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))

  plot(x=RT_df$true,y=RT_df$detected,ylim=c(0,80),xlim=c(0,80),pch=20, xlab="manual RT", ylab="picked RT", main=paste0("Manual annotation benchmark\nTPR = ",round(TPR,digits=3),"\nFPR = ",round(FPR,digits=3)))
  abline(mod1,col="red")
  abline(0,1,lty=2)
  text(x = 19, y = 2.5, labels = mylabel,col="red")
}

estimate_TPR_and_FPR <- function(detected.complexes,true.complexes,input.complexes){
  TP <- sum(detected.complexes %in% true.complexes)
  FN <- length(true.complexes) - TP
  FP <- sum(!(detected.complexes %in% true.complexes))
  negative.complexes <- input.complexes[!(input.complexes %in% true.complexes)]
  TN <- sum(!(negative.complexes %in% detected.complexes))
  
  TPR <- TP / (TP + FN)
  FPR <- FP / (TN + FP)
  list(TPR=TPR,FPR=FPR)
}
