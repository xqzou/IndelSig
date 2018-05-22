## A simple parser for defining and searching for mutation signatures in a processed data-set
## Grammar
# [(+/-)ATCG]
# [(+/-) (ATCG | (>,<,=) 1)](Rep (>,<,=) 1) | (Mh (>,<,=) 1)

## Libraries
library(dplyr)

## Dictionaries defining components 
type_operators <- c('Ins'='+','Del'='-','Complex'='Com')
size_operators <- c('>','<','=')
keywords <- c('Microhomology-mediated'='Mh','Repeat-mediated'='Rep','Ins'='Ins','None'='Others')
####################################################################################################
extract_regions <- function(channel){
  prime5_channel <- regexpr('^.+?\\[',channel)
  indel_channel <- regexpr('\\[(.+?)\\]',channel)
  prime3_channel <- regexpr('(\\].+?$)',channel)
  
  prime5 <- substr(channel,prime5_channel[1],attr(prime5_channel,'match.length')-1)
  indel <- substr(channel,indel_channel[1]+1,attr(indel_channel,'match.length')+indel_channel[1]-2)
  prime3 <- substr(channel,prime3_channel[1]+1,attr(prime3_channel,'match.length')+prime3_channel[1])
  
  result <- list(prime5=prime5,indel=indel,prime3=prime3)
  result <- lapply(result,
                   function(x){
                     y <- strsplit(x,'')[[1]]
                     y <- y[y!=' ']})
  return(result)
}

## Looking up using 'change' rather than 'change_pyr'
indel_parse <- function(str,type_operators,size_operators){
  ## Determine mutation type
  if(str[1] %in% type_operators){
    mutation_type <- names(type_operators)[str[1]==type_operators]
    str <- str[-1]
  }
  else {
    mutation_type <- names(type_operators)
  }

  ## Parse the string to generate queries.
  ## Either a length query, or an indentity query
  query <- ''
  for(i in seq_along(str)){ ## Size query
    if(str[i] %in% size_operators){
      query <- paste('df["indel_length"]')
      for(j in i:length(str)){
        query <- paste0(query,str[j])
      }
      ops <- str[i:(i+1)]
      if('=' %in% ops == T & '>' %in% ops == F & '<' %in% ops == F ){
        query <- gsub(pattern = '=',replacement = '==',x=query)
      }
      break
    }
    else{ ## character string
      if(str[i] %in% c('b','p')){
        next()
      }
      else{
        ## Add in pyrimidine / purine lookup
        query <- 'df["change_pyr"]==\"'
        chars <- ''
        for(j in i:length(str)){
          chars <- paste0(chars,str[j]) 
        }
        query <- paste0(query,chars,'\"')
        break
      }
    }
  }
  return(list(mutation_type = mutation_type,query=query))
}

prime_parse <- function(str,size_operators,keywords,region){
  ## If there's no operators, return an empty query
  if(length(str) == 0){
    query <- ''
  } 
  else {
    ## Look first for the for size operator 
    size_op_check <- unlist(lapply(size_operators,grepl,str))
    
    if(any(size_op_check)){ ## if there is one, set the keyword as the text the left of it
      size_op <- min(which(size_op_check) %% length(str))
      keyword <- paste0(str[1:(size_op-1)],collapse='')
    }
    else{ ## otherwise, set the keyword as the whole string
      keyword <- keyword <- paste0(str,collapse='')
    }
     
    ## If the keyword appears in the dictionary of keywords i.e Mh, Rep or Ins
    if(keyword %in% keywords){
      ## Obtain its searchable name in the data.frame
      keyword <- names(keywords)[keyword == keywords]
      
      ## If there are any size operators next to it, i.e > or <, then
      ## add an additional query for 'repcount' *will be bp in next version*
      if(any(size_op_check)){
        query_1 <- paste0('df[["Classification"]] == "',keyword,'" & ')
        query_2 <- paste0(str[size_op:length(str)],collapse = '')
        ops <- str[size_op:(size_op+1)]
        if('=' %in% ops == T & '>' %in% ops == F & '<' %in% ops == F ){
          query_2 <- gsub(pattern = '=',replacement = '==',x=query_2)
        }
        query <- paste0(query_1,'as.numeric(as.character(df[["Score"]]))',query_2)
      } 
      else{ ## Otherwise, generate a simple query for the Classification
        query <- paste0('df[["Classification"]] == \"',keyword,'\"')
      }
    } 
    else{ ## Otherwise, match the string exactly 
      chars_joined <- paste0(str,collapse = "")
      chars_joined <- paste0('\"',chars_joined,'\"')
      if(region==3){
        query <- paste0('substr(df[["slice',region,'_pyr\"]],1,',length(str),') == ',chars_joined)
      }
      else{
        query <- paste0('substr(df[["slice',region,'_pyr\"]], nchar(df[["slice',region,'_pyr\"]]) - ',length(str),'+1,','nchar(df[["slice',region,'_pyr\"]])) == ',chars_joined)
      }
    }
  }
  result <- list(query= query)
  return(result)
}

parse_channel <- function(df,str,subset=F,id=F,verbose=F){
  #browser()
  ## First generate the regions
  regions <- extract_regions(str)
  
  ## Parse the results into commands
  indel_parsed <- indel_parse(regions$indel,type_operators,size_operators)
  prime3_parsed <- prime_parse(regions$prime3,size_operators,keywords,3)
  prime5_parsed <- prime_parse(regions$prime5,size_operators,keywords,5)
  
  ## Filter the data-set by which mutation type it is. 
  mut_type <- df$Type %in% indel_parsed$mutation_type
  queries <- c(indel_parsed$query,prime3_parsed$query,prime5_parsed$query)
  queries <- queries[queries!='']
  total_query <- paste(queries,collapse = ' & ')
  
  if(verbose) print(c(total_query,indel_parsed$mutation_type))
  
  if(subset){
    df <- df[eval(parse(text=total_query)) & mut_type,]
  }
  else if(id){
    df <- which(eval(parse(text=total_query)) & mut_type)
  }
  else{
    df[eval(parse(text=total_query)) & mut_type,'Subtype'] <- str
  }
  return(df)
}

parse_channels <- function(dataset,channels,subset=F,na.rm=F){
  result <- list()
  for(chan in channels){
    if(subset){
      result[[chan]] <- parse_channel(dataset,chan,subset = subset)
    }
    else{
      dataset <- parse_channel(dataset,chan)
    }
  }
  
  if(na.rm){
    dataset <- dataset[!is.na(dataset$Subtype),]
  }
  
  if(subset){
    return(result)
  } else {
    return(dataset)}
}


assign_mutation_type <- function(df,groups,map=F){
  if(all(grepl(pattern = '\\[*\\]',groups))==F){
    groups <- gsub(pattern = '^',replacement = '[',groups)
    groups <- gsub(pattern = '$',replacement = ']',groups)
  }
  channels <- unique(df$Subtype)
  result <- NULL
  for( i in seq_along(groups)){
    if(map){
      channels[match(df$Subtype[parse_channel(df,groups[i],id = T)],channels)] <- groups[i]
      result <- channels
      names(result) <- unique(df$Subtype)
    }
    else{
      df[parse_channel(df,groups[i],id = T),'Mut_Type'] <- groups[i]
      result <- df
    }
  }
  return(result)
}


convert_to_labels <- function(channels){
  fp <- gsub(pattern = 'Ins',replacement = 'Rep',x = channels)
  for(i in seq_along(fp)){
    par <- extract_regions(fp[i])
    if(all(c('R','e','p')==par$prime3[1:3])){
      if(par$prime3[4]=='='){
        repnum <- as.numeric(paste0(par$prime3[5:length(par$prime3)],collapse = ''))
        p <- indel_parse(par$indel,type_operators,size_operators)$query
        if(repnum==0){
          fp[i] <- gsub(pattern = 'Rep=0','NonRep',fp[i])
        }
        else if(grepl('[channel]',p)){
          k <- regexpr(pattern = '[A-Z]+',text =p)
          repstr <- substr(p,k[1],k[1]+attr(k,'match.length')-1)
          repstr <- paste0(rep(repstr,repnum),collapse = '')
          fp[i] <- gsub(paste0('Rep=',repnum),repstr,fp[i])
        }
      }
    }
  }
  return(fp)
}

## Summarize indels
summarize_channel_frequency <- function(df,...){
  args <- list(...)
  if(length(args)==0){
    sum_arg <- 'n()'
    #stop('no variable selected')
  }
  else{
    all_args <- rep('',length(args))
    for(i in seq_along(args)){
      all_args[i] <- paste0(names(args[i]),' %in% c(\"',paste0(unlist(args[i]),collapse = '\",\"'),'\")')
    }
    sum_arg <- paste0('sum(',paste(all_args,collapse = ' & '),')')
  }
  print(sum_arg)
  indel_catalog <- df %>%
    group_by(Subtype,Mut_Type) %>%
    summarize_('aggregate'=sum_arg) %>%
    .[match(exp_channels,.$Subtype),]
  return(indel_catalog)
}
