library(ggplot2)
library(plyr) 
library(reshape2)
library(gridExtra)

plotbasis_aggragated_indel_29types_6 <- function(muts_basis,h,w,outputname){
  mean_parentmuts <- sum(muts_basis[,1:(dim(muts_basis)[2]-1)])/(dim(muts_basis)[2]-1)
  muts_basis$aggragate <- rowSums(muts_basis[,1:(dim(muts_basis)[2]-1)])/sum(muts_basis[,1:(dim(muts_basis)[2]-1)])
  
  muts_template <- read.table("/nfs/cancer_archive04/xz3/h_mutagen/27_indel_29channels/indel_29channels_6.txt",sep = "\t",header = T, as.is = T)
  muts_basis <- merge(muts_template, muts_basis,by="indelsubtype",all.x=T)
  muts_basis[is.na(muts_basis)] <- 0
  
  # E69F00(light orange, +C); D55E00(dark orange, +T); F0E442(yellow, +>1); 56B4E9(light blue, -C); 0072B2(dark blue, -T); 009E73(green, ->1); CC79A7(pink, MhDel)
  #indel_mypalette <- c("#009E73","#56B4E9","#0072B2","#F0E442","#E69F00","#D55E00")
  
  # 009E73 green, ->1; 56B4E9 light blue -C; E69F00 light orange -T; CC79A7 pink +>1; 0072B2 dark blue +C; D55E00 dark orange +T
  indel_mypalette <- c("#009E73","#56B4E9","#E69F00","#CC79A7","#0072B2","#D55E00")
  
  #[+C]A; [+C]G; [+C]T; [+C]C; [+C]CC; [+C]LongRep (Rep_count>2)
  #[+T]A; [+T]C; [+T]G; [+T]T; [+T]TT; [+T]LongRep (Rep_count>2)
  #[+>1bp]NonRep; [+>1bp]Rep (Rep_count>0)
  #[-C]A; [-C]G; [-C]T; [-C]C; [-C]CC; [-C]LongRep (Rep_count>2)
  #[-T]A; [-T]C; [-T]G; [-T]T; [-T]TT; [-T]LongRep (Rep_count>2)
  #[->1bp]Others; [->1bp]Rep (Rep_count>0); [->1bp]Mh
  
  
  indel_positions <- c("[+C]A","[+C]G","[+C]T","[+C]C","[+C]CC","[+C]LongRep",
                       "[+T]A","[+T]C","[+T]G","[+T]T","[+T]TT","[+T]LongRep",
                       "[+>1]NonRep", "[+>1]Rep",
                       "[-C]A","[-C]G","[-C]T","[-C]C","[-C]CC","[-C]LongRep",
                       "[-T]A","[-T]C","[-T]G","[-T]T","[-T]TT","[-T]LongRep",
                       "[->1]Others","[->1]Rep","[->1]Mh") 
  
  indel_labels <- c("[+C]A","[+C]G","[+C]T","[+C]C","[+C]CC","[+C]LR",
                    "[+T]A","[+T]C","[+T]G","[+T]T","[+T]TT","[+T]LR",
                    "[+>1]NonR", "[+>1]Rep",
                    "[-C]A","[-C]G","[-C]T","[-C]C","[-C]CC","[-C]LR",
                    "[-T]A","[-T]C","[-T]G","[-T]T","[-T]TT","[-T]LR",
                    "[->1]Oth.","[->1]Rep","[->1]Mh") 
  
  filename <- paste0(outputname, ".pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis, aes(x=indelsubtype, y=aggragate,fill=indeltype))+ geom_bar(stat="identity",position="dodge")+xlab("Indel Types")+ylab("Percentage")
  p <- p+scale_y_continuous(limits=c(0,1),breaks=(seq(0,1,0.2)),labels = scales::percent)
  p <- p+scale_x_discrete(limits = indel_positions,labels = indel_labels)+ggtitle(outputname)
  p <- p+scale_fill_manual(values=indel_mypalette)
  p <- p+theme(axis.text.x=element_text(angle=45, vjust=0.5, size=10,colour = "black"),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
  
  print(p)
  dev.off()
}
plotCountbasis_aggragated_indel_29types_6 <- function(muts_basis,h,w,outputname){
  mean_parentmuts <- sum(muts_basis[,1:(dim(muts_basis)[2]-1)])/(dim(muts_basis)[2]-1)
  muts_basis$aggragate <- rowSums(muts_basis[,1:(dim(muts_basis)[2]-1)])
  
  muts_template <- read.table("/nfs/cancer_archive04/xz3/h_mutagen/27_indel_29channels/indel_29channels_6.txt",sep = "\t",header = T, as.is = T)
  muts_basis <- merge(muts_template, muts_basis,by="indelsubtype",all.x=T)
  muts_basis[is.na(muts_basis)] <- 0
  
  # E69F00(light orange, +C); D55E00(dark orange, +T); F0E442(yellow, +>1); 56B4E9(light blue, -C); 0072B2(dark blue, -T); 009E73(green, ->1); CC79A7(pink, MhDel)
  #indel_mypalette <- c("#009E73","#56B4E9","#0072B2","#F0E442","#E69F00","#D55E00")
  
  # 009E73 green, ->1; 56B4E9 light blue -C; E69F00 light orange -T; CC79A7 pink +>1; 0072B2 dark blue +C; D55E00 dark orange +T
  indel_mypalette <- c("#009E73","#56B4E9","#E69F00","#CC79A7","#0072B2","#D55E00")
  
  #[+C]A; [+C]G; [+C]T; [+C]C; [+C]CC; [+C]LongRep (Rep_count>2)
  #[+T]A; [+T]C; [+T]G; [+T]T; [+T]TT; [+T]LongRep (Rep_count>2)
  #[+>1bp]NonRep; [+>1bp]Rep (Rep_count>0)
  #[-C]A; [-C]G; [-C]T; [-C]C; [-C]CC; [-C]LongRep (Rep_count>2)
  #[-T]A; [-T]C; [-T]G; [-T]T; [-T]TT; [-T]LongRep (Rep_count>2)
  #[->1bp]Others; [->1bp]Rep (Rep_count>0); [->1bp]Mh
  
  
  indel_positions <- c("[+C]A","[+C]G","[+C]T","[+C]C","[+C]CC","[+C]LongRep",
                       "[+T]A","[+T]C","[+T]G","[+T]T","[+T]TT","[+T]LongRep",
                       "[+>1]NonRep", "[+>1]Rep",
                       "[-C]A","[-C]G","[-C]T","[-C]C","[-C]CC","[-C]LongRep",
                       "[-T]A","[-T]C","[-T]G","[-T]T","[-T]TT","[-T]LongRep",
                       "[->1]Others","[->1]Rep","[->1]Mh") 
  
  indel_labels <- c("[+C]A","[+C]G","[+C]T","[+C]C","[+C]CC","[+C]LR",
                    "[+T]A","[+T]C","[+T]G","[+T]T","[+T]TT","[+T]LR",
                    "[+>1]NonR", "[+>1]Rep",
                    "[-C]A","[-C]G","[-C]T","[-C]C","[-C]CC","[-C]LR",
                    "[-T]A","[-T]C","[-T]G","[-T]T","[-T]TT","[-T]LR",
                    "[->1]Oth.","[->1]Rep","[->1]Mh") 
  
  filename <- paste0(outputname, ".pdf")
  pdf(file=filename, onefile=TRUE,width=w,height=h)
  p <- ggplot(data=muts_basis, aes(x=indelsubtype, y=aggragate,fill=indeltype))+ geom_bar(stat="identity",position="dodge")+xlab("Indel Types")+ylab("Percentage")
  #  p <- p+scale_y_continuous(limits=c(0,40),breaks=(seq(0,40,10)))
  p <- p+scale_x_discrete(limits = indel_positions,labels = indel_labels)+ggtitle(outputname)
  p <- p+scale_fill_manual(values=indel_mypalette)
  p <- p+theme(axis.text.x=element_text(angle=45, vjust=0.5, size=10,colour = "black"),
               axis.text.y=element_text(size=10,colour = "black"),
               axis.title.x = element_text(size=15),
               axis.title.y = element_text(size=15),
               plot.title = element_text(size=10),
               panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank(),
               panel.grid.major.y = element_blank(),
               panel.grid.minor.y = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.border = element_rect(colour = "black", fill=NA))
  
  print(p)
  dev.off()
  
  return(muts_basis)
}

samples_details <- read.table("/nfs/cancer_archive04/xz3/h_mutagen/00_common/mutagen_info_forR.txt",sep = "\t",header = T,as.is = T,quote = "\"")
indel.classified <- read.table("./indel.classified.txt",sep = "\t",header = T, as.is = T)
indel.classified$Sample.Name <- sub("\\_.*","",indel.classified$Sample)
indel.classified <- indel.classified[indel.classified$Sample.Name!="MSM0",]

indel.classified_details <- merge(indel.classified, samples_details, by="Sample.Name")
indel.classified_details <- indel.classified_details[indel.classified_details$VAF.Tum_Cal>=0.2,]


###############################################
#
# Indel catalogue 29 channels （ 27_indel_29channels）

#[+C]A; [+C]G; [+C]T; [+C]C; [+C]CC; [+C]LongRep (Rep_count>2)
#[+T]A; [+T]C; [+T]G; [+T]T; [+T]TT; [+T]LongRep (Rep_count>2)
#[+>1bp]NonRep; [+>1bp]Rep (Rep_count>0)
#[-C]A; [-C]G; [-C]T; [-C]C; [-C]CC; [-C]LongRep (Rep_count>2)
#[-T]A; [-T]C; [-T]G; [-T]T; [-T]TT; [-T]LongRep (Rep_count>2)
#[->1bp]Others; [->1bp]Rep (Rep_count>0); [->1bp]Mh

###############################################
indel.classified_details2 <- indel.classified_details
indel.classified_details2$Subtype <- NULL
indel.classified_details2[indel.classified_details2$Type=="Ins" & indel.classified_details2$change.pyr=="C" & indel.classified_details2$repcount==0 & indel.classified_details2$slice3_1bp_pyr=="A","Subtype"] <- "[+C]A"
indel.classified_details2[indel.classified_details2$Type=="Ins" & indel.classified_details2$change.pyr=="C" & indel.classified_details2$repcount==0 & indel.classified_details2$slice3_1bp_pyr=="G","Subtype"] <- "[+C]G"
indel.classified_details2[indel.classified_details2$Type=="Ins" & indel.classified_details2$change.pyr=="C" & indel.classified_details2$repcount==0 & indel.classified_details2$slice3_1bp_pyr=="T","Subtype"] <- "[+C]T"
indel.classified_details2[indel.classified_details2$Type=="Ins" & indel.classified_details2$change.pyr=="C" & indel.classified_details2$repcount==1,"Subtype"] <- "[+C]C"
indel.classified_details2[indel.classified_details2$Type=="Ins" & indel.classified_details2$change.pyr=="C" & indel.classified_details2$repcount==2,"Subtype"] <- "[+C]CC"
indel.classified_details2[indel.classified_details2$Type=="Ins" & indel.classified_details2$change.pyr=="C" & indel.classified_details2$repcount>2,"Subtype"] <- "[+C]LongRep"

indel.classified_details2[indel.classified_details2$Type=="Ins" & indel.classified_details2$change.pyr=="T" & indel.classified_details2$repcount==0 & indel.classified_details2$slice3_1bp_pyr=="A","Subtype"] <- "[+T]A"
indel.classified_details2[indel.classified_details2$Type=="Ins" & indel.classified_details2$change.pyr=="T" & indel.classified_details2$repcount==0 & indel.classified_details2$slice3_1bp_pyr=="C","Subtype"] <- "[+T]C"
indel.classified_details2[indel.classified_details2$Type=="Ins" & indel.classified_details2$change.pyr=="T" & indel.classified_details2$repcount==0 & indel.classified_details2$slice3_1bp_pyr=="G","Subtype"] <- "[+T]G"
indel.classified_details2[indel.classified_details2$Type=="Ins" & indel.classified_details2$change.pyr=="T" & indel.classified_details2$repcount==1,"Subtype"] <- "[+T]T"
indel.classified_details2[indel.classified_details2$Type=="Ins" & indel.classified_details2$change.pyr=="T" & indel.classified_details2$repcount==2,"Subtype"] <- "[+T]TT"
indel.classified_details2[indel.classified_details2$Type=="Ins" & indel.classified_details2$change.pyr=="T" & indel.classified_details2$repcount>2,"Subtype"] <- "[+T]LongRep"

indel.classified_details2[indel.classified_details2$Type=="Ins" & indel.classified_details2$indel.length>1 & indel.classified_details2$repcount<1,"Subtype"] <- "[+>1]NonRep"
indel.classified_details2[indel.classified_details2$Type=="Ins" & indel.classified_details2$indel.length>1 & indel.classified_details2$repcount>=1,"Subtype"] <- "[+>1]Rep"

indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$change.pyr=="C" & indel.classified_details2$repcount==0 & indel.classified_details2$slice3_1bp_pyr=="A","Subtype"] <- "[-C]A"
indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$change.pyr=="C" & indel.classified_details2$repcount==0 & indel.classified_details2$slice3_1bp_pyr=="G","Subtype"] <- "[-C]G"
indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$change.pyr=="C" & indel.classified_details2$repcount==0 & indel.classified_details2$slice3_1bp_pyr=="T","Subtype"] <- "[-C]T"
indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$change.pyr=="C" & indel.classified_details2$repcount==1,"Subtype"] <- "[-C]C"
indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$change.pyr=="C" & indel.classified_details2$repcount==2,"Subtype"] <- "[-C]CC"
indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$change.pyr=="C" & indel.classified_details2$repcount>2,"Subtype"] <- "[-C]LongRep"


#indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$classification=="Repeat-mediated" & indel.classified_details2$change.pyr=="C" & indel.classified_details2$repcount>1,"Subtype"] <- "[-C]_Rep"

indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$change.pyr=="T" & indel.classified_details2$repcount==0 & indel.classified_details2$slice3_1bp_pyr=="A","Subtype"] <- "[-T]A"
indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$change.pyr=="T" & indel.classified_details2$repcount==0 & indel.classified_details2$slice3_1bp_pyr=="C","Subtype"] <- "[-T]C"
indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$change.pyr=="T" & indel.classified_details2$repcount==0 & indel.classified_details2$slice3_1bp_pyr=="G","Subtype"] <- "[-T]G"
indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$change.pyr=="T" & indel.classified_details2$repcount==1,"Subtype"] <- "[-T]T"
indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$change.pyr=="T" & indel.classified_details2$repcount==2,"Subtype"] <- "[-T]TT"
indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$change.pyr=="T" & indel.classified_details2$repcount>2,"Subtype"] <- "[-T]LongRep"
#indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$classification=="Repeat-mediated" & indel.classified_details2$change.pyr=="T" & indel.classified_details2$repcount>1,"Subtype"] <- "[-T]_Rep"

indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$classification !="Microhomology-mediated" & indel.classified_details2$indel.length>1 & indel.classified_details2$repcount<1,"Subtype"] <- "[->1]Others"
indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$classification=="Repeat-mediated" & indel.classified_details2$indel.length>1 & indel.classified_details2$repcount>=1,"Subtype"] <- "[->1]Rep"
indel.classified_details2[indel.classified_details2$Type=="Del" & indel.classified_details2$classification =="Microhomology-mediated","Subtype"] <- "[->1]Mh"



indel_catalogue <- data.frame(table(indel.classified_details2$Sample,indel.classified_details2$Subtype))
names(indel_catalogue) <- c("suclone","subtype","freq")
indel_catalogue <- dcast(indel_catalogue,suclone~subtype)
indel_catalogue$Sample.Name <- sub("\\_.*","",indel_catalogue$suclone)
write.table(indel_catalogue, paste0("indel_catalogue",".txt"),sep = "\t",col.names = T, row.names = F, quote = T)

indel_catalogue <- merge(indel_catalogue,samples_details, by="Sample.Name", all.x=T)
indel_catalogue$indel_num <- rowSums(indel_catalogue[,3:(2+29)])
write.table(indel_catalogue, paste0("indel_catalogue_withinfo",".txt"),sep = "\t",col.names = T, row.names = F, quote = F)

muts_control <- indel_catalogue[indel_catalogue$Group=="Control",]
muts_compound <- indel_catalogue[!indel_catalogue$Group=="Control",]

###############################################
# Plot indel profiles
###############################################
if(TRUE){
  indel_catalogue <- read.table("indel_catalogue_withinfo.txt",sep = "\t",header = T, as.is = T, quote = "\"",check.names = F)
  muts_control <- indel_catalogue[indel_catalogue$Group=="Control",]
  muts_compound <- indel_catalogue[indel_catalogue$Group!="Control",]
  
  # Control
  parentmuts <- t(muts_control[,3:(2+29)])
  colnames(parentmuts) <- muts_control$suclone
  parentmuts <- as.data.frame(parentmuts)
  
  parentmuts$indelsubtype <- rownames(parentmuts)
  control_mean <- plotCountbasis_aggragated_indel_29types_6(parentmuts,3,10,"controlbasis_indels_Count")
  plotbasis_aggragated_indel_29types_6(parentmuts,3,10,"controlbasis_indels")
  
  # All compound
  compoundlist <- data.frame(table(muts_compound$Sample.Name))
  names(compoundlist) <- c("Sample.Name","Freq")
  child_mean_all <- NULL
  for(i in 1:dim(compoundlist)[1]){
    
    print(i)
    currentcompound <- muts_compound[muts_compound$Sample.Name==as.character(compoundlist[i,1]),]
    
    currentcompound_muts <- currentcompound[,3:(2+29)]
    currentcompound_muts <- t(currentcompound_muts)
    colnames(currentcompound_muts) <- currentcompound$suclone
    currentcompound_muts <- as.data.frame(currentcompound_muts)
    currentcompound_muts$indelsubtype <- rownames(currentcompound_muts)
    currentcompound_mean <- plotCountbasis_aggragated_indel_29types_6(currentcompound_muts,3,10,paste0("indels_",compoundlist[i,1],"_",currentcompound[1,"Compound.Abbreviation"],"_Count"))
    write.table(currentcompound_mean,paste0("indels_",compoundlist[i,1],"_",currentcompound[1,"Compound.Abbreviation"],"_count.txt"),sep = "\t",col.names = T, row.names = F, quote = F)
    plotbasis_aggragated_indel_29types_6(currentcompound_muts,3,10,paste0("indels_",compoundlist[i,1],"_",currentcompound[1,"Compound.Abbreviation"],"_Percentage"))
    
    child_mean_all <- cbind(child_mean_all,currentcompound_mean$aggragate)
  }
  
  child_mean_all <- data.frame(child_mean_all)
  names(child_mean_all) <- compoundlist$Sample.Name
  
  control_child_mean <- cbind(control_mean$aggragate,child_mean_all)
  control_child_mean$indelsubtype <- control_mean$indelsubtype
  write.table(control_child_mean,"control_child_aggragate.txt",sep = "\t",col.names = T, row.names = F, quote = F)
}


