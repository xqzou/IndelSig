################################################################################
## Main function
################################################################################
prepare_indel_data <- function(df,ref='Ref',alt='Alt',type='Type',chrom='Chrom',pos='Pos',mut_type= c(Ins='Ins',Del='Del',Complex='Complex')) {
  args <- as.list(environment(),all=TRUE)[-c(1,7)]
  if (nrow(df)<1 | all(unlist(args) %in% colnames(df)) != T) {
    stop('Dataframe is invalid. Please check column names.')
  }
  
  ## Determine mutation change
  df$change <- NULL
  mut_ref <- c('Alt','Ref','Ref')
  for(i in seq_along(mut_type)){
    cur_mut <- df[[type]] == mut_type[[i]]
    df[cur_mut,'change'] <- substr(as.character(df[cur_mut, mut_ref[i]]),2,1e5)
  }
  df$indel_length <- nchar(df$change)
  
  ## Determine 5' and 3' contexts.
  ## length of indel is included because of long repeat insertion / deletions exceding exactly 30 bp upstream / downstream context. 
  ## Determine the 5' context
  df$extend5 <- df[[pos]]- df$indel_length -30
  df$slice5 <- as.character(getSeq(Hsapiens, paste0('chr',df[[chrom]]), df$extend5, df[[pos]]))
  
  ## Determine the 3' Context
  df$extend3 = df[[pos]] + 2*df$indel_length + 30 # At least 30 bp upstream, plus adjustment for deletion events
  ## For Deletions
  df[df[[type]]==mut_type['Del'],"slice3"] <- as.character(getSeq(Hsapiens,
                                                               paste0('chr',df[df[[type]]==mut_type['Del'],chrom]),
                                                               df[df[[type]]==mut_type['Del'],pos] + df[df[[type]]==mut_type['Del'],'indel_length'] + 1,
                                                               df[df[[type]]==mut_type['Del'],'extend3']))
  ## For insertions
  df[df[[type]]==mut_type['Ins'],"slice3"] <- as.character(getSeq(Hsapiens, 
                                                               paste0('chr', df[df[[type]]==mut_type['Ins'],chrom]),
                                                               df[df[[type]]==mut_type['Ins'],pos]+1,
                                                               df[df[[type]]==mut_type['Ins'],"extend3"]))
  ## Pyrimide base changes. 
  ## Single base changes
  df$change_pyr <- df$change
  df$slice3_pyr <- df$slice3
  df$slice5_pyr <- df$slice5
  
  pyr_bases <- c('T'='A','C'='G')
  
  for(i in seq_along(pyr_bases)){
    cur_base <- df$change == pyr_bases[i]
    df[cur_base,'change_pyr'] <- names(pyr_bases)[i]
    df[cur_base,'slice5_pyr'] <- as.character(reverseComplement(DNAStringSet(df[cur_base,'slice3'])))
    df[cur_base,'slice3_pyr'] <- as.character(reverseComplement(DNAStringSet(df[cur_base,'slice5'])))
  }
  
  return(df)
}