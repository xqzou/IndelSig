# Dominik
# Modified by Xueqing

# prepare indel in tab file,
prepare.indel.df.tab <- function(indel.data) {
  
  if (nrow(indel.data)>0) {
    
    
    indel.data$ref.length <- nchar(indel.data$Ref)
    indel.data$alt.length <- nchar(indel.data$Alt)
    indel.data$indel.length <- abs(indel.data$ref.length - indel.data$alt.length)
    
     
    # sequence of change
    indel.data$change <- NULL
    indel.data[indel.data$Type=="Complex","change"] <- substr( as.character(indel.data[indel.data$Type=='Complex',"Ref"]),2,1e5)
    indel.data[indel.data$Type=="Ins","change"] <- substr( as.character(indel.data[indel.data$Type=='Ins',"Alt"]),2,1e5)
    indel.data[indel.data$Type=="Del","change"] <- substr( as.character(indel.data[indel.data$Type=='Del',"Ref"]),2,1e5)
     
    
    # 27bp before and after change                     
    indel.data$extend5 = indel.data$Pos-indel.data$indel.length-25;
    indel.data$extend3 = indel.data$Pos+indel.data$indel.length + indel.data$indel.length+25;
    
    
    indel.data$slice5 <- as.character(getSeq(Hsapiens, paste0('chr',indel.data$Chrom), indel.data$extend5, indel.data$Pos))
    indel.data$slice3 <- NULL
    indel.data[indel.data$Type=="Del","slice3"] <- as.character(getSeq(Hsapiens, paste0('chr',indel.data[indel.data$Type=="Del","Chrom"]), indel.data[indel.data$Type=="Del","Pos"]+indel.data[indel.data$Type=="Del","indel.length"]+1, indel.data[indel.data$Type=="Del","extend3"]))
    indel.data[indel.data$Type=="Ins","slice3"] <- as.character(getSeq(Hsapiens, paste0('chr',indel.data[indel.data$Type=="Ins","Chrom"]), indel.data[indel.data$Type=="Ins","Pos"]+1, indel.data[indel.data$Type=="Ins","extend3"]))
    
    # 1bp before and after change                     
    indel.data$slice5_1bp <- as.character(getSeq(Hsapiens, paste0('chr',indel.data$Chrom), indel.data$Pos, indel.data$Pos))
    indel.data$slice3_1bp <- NULL
    indel.data[indel.data$Type=="Del","slice3_1bp"] <- as.character(getSeq(Hsapiens, paste0('chr',indel.data[indel.data$Type=="Del","Chrom"]), indel.data[indel.data$Type=="Del","Pos"]+indel.data[indel.data$Type=="Del","indel.length"]+1, indel.data[indel.data$Type=="Del","Pos"]+indel.data[indel.data$Type=="Del","indel.length"]+1))
    indel.data[indel.data$Type=="Ins","slice3_1bp"] <- as.character(getSeq(Hsapiens, paste0('chr',indel.data[indel.data$Type=="Ins","Chrom"]), indel.data[indel.data$Type=="Ins","Pos"]+1, indel.data[indel.data$Type=="Ins","Pos"]+1))
    
    # 2bp before and after change                     
    
    indel.data$slice5_2bp <- as.character(getSeq(Hsapiens, paste0('chr',indel.data$Chrom), indel.data$Pos-1, indel.data$Pos))
    indel.data$slice3_2bp <- NULL
    indel.data[indel.data$Type=="Del","slice3_2bp"] <- as.character(getSeq(Hsapiens, paste0('chr',indel.data[indel.data$Type=="Del","Chrom"]), indel.data[indel.data$Type=="Del","Pos"]+indel.data[indel.data$Type=="Del","indel.length"]+1, indel.data[indel.data$Type=="Del","Pos"]+indel.data[indel.data$Type=="Del","indel.length"]+2))
    indel.data[indel.data$Type=="Ins","slice3_2bp"] <- as.character(getSeq(Hsapiens, paste0('chr',indel.data[indel.data$Type=="Ins","Chrom"]), indel.data[indel.data$Type=="Ins","Pos"]+1, indel.data[indel.data$Type=="Ins","Pos"]+2))
    
    
    # 3bp before and after change                     
      
#    indel.data$slice5_3bp <- as.character(getSeq(Hsapiens, paste0('chr',indel.data$Chrom), indel.data$Pos-2, indel.data$Pos))
#    indel.data$slice3_3bp <- NULL
#    indel.data[indel.data$Type=="Del","slice3_3bp"] <- as.character(getSeq(Hsapiens, paste0('chr',indel.data[indel.data$Type=="Del","Chrom"]), indel.data[indel.data$Type=="Del","Pos"]+indel.data[indel.data$Type=="Del","indel.length"]+1, indel.data[indel.data$Type=="Del","Pos"]+indel.data[indel.data$Type=="Del","indel.length"]+3))
#    indel.data[indel.data$Type=="Ins","slice3_3bp"] <- as.character(getSeq(Hsapiens, paste0('chr',indel.data[indel.data$Type=="Ins","Chrom"]), indel.data[indel.data$Type=="Ins","Pos"]+1, indel.data[indel.data$Type=="Ins","Pos"]+3))
    
    # Pyrimidine represnetation for 1bp indels
    indel.data$change.pyr <- indel.data$change
    indel.data$slice3_1bp_pyr <- indel.data$slice3_1bp
    indel.data$slice5_1bp_pyr <- indel.data$slice5_1bp
    
    indel.data$slice3_2bp_pyr <- indel.data$slice3_2bp
    indel.data$slice5_2bp_pyr <- indel.data$slice5_2bp
    
    indel.data[indel.data$change=="A","change.pyr"] <- "T"
    indel.data[indel.data$change=="A","slice5_1bp_pyr"] <- as.character(reverseComplement(DNAStringSet(indel.data[indel.data$change=="A","slice3_1bp"])))
    indel.data[indel.data$change=="A","slice3_1bp_pyr"] <- as.character(reverseComplement(DNAStringSet(indel.data[indel.data$change=="A","slice5_1bp"])))
    
    indel.data[indel.data$change=="A","slice5_2bp_pyr"] <- as.character(reverseComplement(DNAStringSet(indel.data[indel.data$change=="A","slice3_2bp"])))
    indel.data[indel.data$change=="A","slice3_2bp_pyr"] <- as.character(reverseComplement(DNAStringSet(indel.data[indel.data$change=="A","slice5_2bp"])))
    
    
    indel.data[indel.data$change=="G","change.pyr"] <- "C"
    indel.data[indel.data$change=="G","slice5_1bp_pyr"] <- as.character(reverseComplement(DNAStringSet(indel.data[indel.data$change=="G","slice3_1bp"])))
    indel.data[indel.data$change=="G","slice3_1bp_pyr"] <- as.character(reverseComplement(DNAStringSet(indel.data[indel.data$change=="G","slice5_1bp"])))
    
    indel.data[indel.data$change=="G","slice5_2bp_pyr"] <- as.character(reverseComplement(DNAStringSet(indel.data[indel.data$change=="G","slice3_2bp"])))
    indel.data[indel.data$change=="G","slice3_2bp_pyr"] <- as.character(reverseComplement(DNAStringSet(indel.data[indel.data$change=="G","slice5_2bp"])))
    
    
    # indel.df needs following columns:
    # indel.type
    # change
    # slice3
    # slice5
    # indel.length
    
    return(indel.data)
    

  }
}
