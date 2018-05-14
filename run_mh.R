# Dominik Glodzik
source("https://bioconductor.org/biocLite.R")
biocLite("VariantAnnotation")
biocLite("GenomicFeatures")
biocLite("BSgenome.Hsapiens.UCSC.hg19")

library("VariantAnnotation")
library('GenomicFeatures')
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)

source('./microhomology.R')
source('./prepare.indel.df.R')


indel.data <- read.table("/nfs/cancer_archive04/xz3/h_mutagen/22_indels/00_data/denovo_subclone_indels_final_freeze.txt",sep="\t",header = T, as.is = T)
# No complex
indel.data <- indel.data[indel.data$Type %in% c("Del","Ins"),]

# convert formats, and find context of the indels
indel.df <- prepare.indel.df.tab(indel.data)
indel.df.max100 <- indel.df[indel.df$indel.length<=100,]
# indel classification
indel.classified.df <- mh_indel(indel.df.max100)
write.table(indel.classified.df, "indel.classified.txt",sep = "\t",col.names = T, row.names = F, quote = F)


