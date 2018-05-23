## Sample workflow for processing an Indel data-set
## Load Bioconductor + Repos
# source("https://bioconductor.org/biocLite.R")
# biocLite("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)
library(magrittr)

## Load the caller and helper functions
source('./v2/prepare_indel_dataframe.R')
source('./v2/callers.R')

## Read in the data-set
indel_data <- read.table("denovo_subclone_indels_final_freeze.txt",sep="\t",header = T, as.is = T)

## Filter to remove Complex mutations
indel_data <- indel_data[indel_data$Type %in% c("Del","Ins"),]

## Prepare the data-set, i.e determine the change, and 5' / 3' sequence contexts
indel_data_prepared <- prepare_indel_data(indel_data)

## Select only 'small' indels.
indel_data_prepared_max100 <- indel_data_prepared[indel_data_prepared$indel_length <= 100,]

## Generate a classified data-set
indel_data_classified <- call_indels(df = indel_data_prepared_max100,
                                     decision = current_decision,
                                      callers = list(microhomology_mediated_deletion,
                                                     repeat_mediated_deletion,
                                                     repeat_insertion))

## in Pipe form!
indel_data <- read.table("denovo_subclone_indels_final_freeze.txt",sep="\t",header = T, as.is = T) %>%
              .[.$Type %in% c('Del','Ins'),] %>%
              prepare_indel_data() %>%
              .[.$indel_length <= 100,] %>%
              call_indels(decision = current_decision,
                          callers = list(microhomology_mediated_deletion,
                                         repeat_mediated_deletion,
                                         repeat_insertion))

## Write out the result
write.table(indel_data_classified, "indel_classified.txt",sep = "\t",col.names = T, row.names = F, quote = F)


