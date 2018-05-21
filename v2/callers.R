################################################################################
## Main function
################################################################################
call_indels <- function(df,decision,callers=list(),top_n = 1){
  channels <- matrix(nrow=nrow(df),ncol = 3*top_n)
  for(i in 1:nrow(df)){ # For each sample
    calls <- c()
    for(caller in callers){ # For each extractor function
      if(df[i,]$Type == attr(caller,'Type')){
        calls <- c(calls, caller(df[i,]) ) # Execute the extractor on the mutation to get a call
      }
    }
    ## Combine all the calls, and execute the decision function to make a final call.
    calls <- matrix(calls,nrow=3)
    colnames(calls) <- calls[1,]
    rownames(calls) <- c('call','score','seq')
    channels[i,1:3] <- decision(calls,top_n) # Evaluate the calls 
  }
  colnames(channels) <- c('Classification','Score','Seq')
  df <- cbind(df,channels)
  return(df)
}

################################################################################
## Decision Functions
################################################################################
current_decision <- function(inputs,top_n=1){
  scores <- inputs['score',]
  mode(scores) <- 'numeric'
  ## If only call is made and it's nonzero, return it.
  if(is.null(inputs)){
    return(cbind(call="None",score=0,seq=''))
  }
  else if(length(scores) == 1 & scores[1] > 0){
    return(inputs)
  }
  ## If working with Del calls
  else if(all(names(scores) %in% c('Microhomology-mediated','Repeat-mediated')) ){
    if(all(scores!=0)){  ## If all scores are nonzero
      ## If more repeat bases than microhomology bases, return Repeat-mediated
      if(scores['Repeat-mediated']* nchar(inputs['seq','Repeat-mediated']) >= scores['Microhomology-mediated']){
        return(inputs[,'Repeat-mediated'])
      }
      else{ ## Else return Mh
        return(inputs[,'Microhomology-mediated'])
      }
    }
    else if (sum(scores) > 0){  ## If one score is nonzero, return that one. 
      return(inputs[,which(scores!=0)])
    }
    else{    ## If both are zero, return a None call
      return(cbind(call="None",score=0,seq=''))
    }
  }
  ## If working with Ins calls return Ins
  else if(all(names(scores) %in% c('Ins'))){
    return(inputs[,'Ins'])
  }
}

################################################################################
## Caller Functions
################################################################################
## Microhomology deletions
microhomology_mediated_deletion <- function(inp){
  ## Use  helper function to find homology between the deletion and 3' context from the left
  mh_seq <- left_shared_subsequence(inp['change'],inp['slice3'])
  result <- c(call="Microhomology-mediated",score=nchar(mh_seq),seq=mh_seq)
  return(result)
}
attr(microhomology_mediated_deletion,'Type') <- 'Del'

## Repat-mediated deletions
repeat_mediated_deletion <- function(inp){
  if(inp['indel_length'] == 1){
    rpt_seq <- left_shared_tandem_subsequence(inp['change'],inp['slice3']) ## should this be +1?
    result <- c(call='Repeat-mediated',
                score=nchar(rpt_seq),
                seq=rpt_seq)
  }
  ## The old style, trying to replicate it without it being ugly 
  else if(inp['change'] == substr(inp['slice3'],1,inp['indel_length'])){
    smallest_rep_unit <- smallest_repetitive_subsequence(inp['change'])
    if (nchar(smallest_rep_unit) < nchar(inp['change'])){
      match <- left_shared_tandem_subsequence(smallest_rep_unit,inp['slice3'])
    }
    else{
      match <- left_shared_tandem_subsequence(inp['change'],inp['slice3'])
    }
    result <- c(call='Repeat-mediated',
                score=nchar(match)/nchar(smallest_rep_unit),
                seq=match)
  }
  else {
    ## Not using old 'chopping' method.
    result <- c(call='Repeat-mediated',
                score=0,
                seq='')
    
  }
  return(result)
}
attr(repeat_mediated_deletion,'Type') <- 'Del'

## Repeat Insertions
repeat_insertion <- function(inp){
  result <- repeat_mediated_deletion(inp)
  result[['call']] <- 'Ins'
  return(result)
}
attr(repeat_insertion,'Type') <- 'Ins'

## Simple test function
dummy_caller <- function(inp){
  return(c(call='Repeat-mediated',score=0,seq='D'))
}
attr(dummy_caller,'Type') <- 'Del'



