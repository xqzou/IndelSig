## Sample workflow for testing out a set of channels and investigating how they vary.
## Load the channel parser and plotting functions
source('v2/parser.R')
source('v2/plot_indel.R')

## First load the data 
samples_details <- read.table("./mutagen_info_forR.txt",sep = "\t",header = T,as.is = T,quote = "\"")
indel_data_classified <- read.table("./indel_classified.txt",sep = "\t",header = T, as.is = T)
indel_data_classified$Sample.Name <- sub("\\_.*","",indel_data_classified$Sample)

## Filter by MSM0
indel_data_classified <- indel_data_classified[indel_data_classified$Sample.Name!="MSM0",]
indel_data_classified <- merge(indel_data_classified, samples_details, by="Sample.Name")

## Filter by VAF.cal > 0.
indel_data_classified <- indel_data_classified[indel_data_classified$VAF.Tum_Cal>=0.2,]



#### Testing a set of channels! #####

## Set up a test set of channels
exp_types <- c('+C','+T','-C','-T','->1','+>1')
exp_channels <- c('[+C]A', '[+C]G','[+C]T', '[+C]Ins=1', '[+C]Ins=2','[+C]Ins>2',
              '[+T]A','[+T]C','[+T]G','[+T]Ins=1','[+T]Ins=2','[+T]Ins>2',
              '[-C]A','[-C]G','[-C]T','[-C]Rep=1','[-C]Rep=2','[-C]Rep>2',
              '[-T]A', '[-T]C', '[-T]G', '[-T]Rep=1', '[-T]Rep=2','[-T]Rep>2',
              '[+>1]Ins=0', '[+>1]Ins>0',
              '[->1]Others', '[->1]Rep>0','[->1]Mh')

## Parse the channels and assign them based on indel classifications 
## then assign a mutation type from 'exp_types', using pipes. 
## Then summarize the data over the control group and plot it
indel_catalog <- indel_data_classified %>% 
                 parse_channels(exp_channels,na.rm = T) %>%
                 assign_mutation_type(exp_types) %>%
                 summarize_channel_frequency(Group='PAH') %>%
                 feature_barplot(x_labels=convert_to_labels(.$Subtype))

