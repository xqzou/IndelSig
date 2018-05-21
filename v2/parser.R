## A simple parser for defining and searching for mutation signatures in a processed data-set
## Grammar
# [(+/-)ATCG]
# [(+/-) (ATCG | (>,<,=) 1)](Rep (>,<,=) 1) | (Mh (>,<,=) 1)

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
      query <- paste('df["indel.length"]')
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
        query <- 'df["change"]==\"'
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
        query_1 <- paste0('df[["classification"]] == "',keyword,'" & ')
        query_2 <- paste0(str[size_op:length(str)],collapse = '')
        ops <- str[size_op:(size_op+1)]
        if('=' %in% ops == T & '>' %in% ops == F & '<' %in% ops == F ){
          query_2 <- gsub(pattern = '=',replacement = '==',x=query_2)
        }
        query <- paste0(query_1,'as.numeric(df[["repcount"]])',query_2)
      } 
      else{ ## Otherwise, generate a simple query for the classification
        query <- paste0('df[["classification"]] == \"',keyword,'\"')
      }
    } 
    else{ ## Otherwise, match the string exactly 
      chars_joined <- paste0(str,collapse = "")
      chars_joined <- paste0('\"',chars_joined,'\"')
      if(region==3){
        query <- paste0('substr(df[["slice',region,'\"]],1,',length(str),') == ',chars_joined)
      }
      else{
        query <- paste0('substr(df[["slice',region,'\"]], nchar(df[["slice',region,'\"]]) - ',length(str),'+1,','nchar(df[["slice',region,'\"]])) == ',chars_joined)
      }
    }
  }
  result <- list(query= query)
  return(result)
}

parse_channel <- function(df,str,show_queries=F,subset=F){
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
  
  if(show_queries) print(total_query)
  
  if(subset){
    df <- df[eval(parse(text=total_query)) & mut_type,]
  }
  else{
    df[eval(parse(text=total_query)) & mut_type,'Subtype'] <- str
  }
  return(df)
}

parse_channels <- function(channels,dataset,subset=F){
  result <- list()
  for(chan in channels){
    if(subset){
      result[[chan]] <- parse_channel(dataset,chan,subset = subset)
    }
    else{
      dataset <- parse_channel(dataset,chan)
    }
  }
  
  if(subset){return(result)} else {return(dataset)}
}





