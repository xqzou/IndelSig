################################################################################
## Helper functions for seq processing
################################################################################
## Determines the shared subsequence between two sequences starting at the left
left_shared_subsequence <- function(seq1,seq2){
  for(i in 1L:min(nchar(seq1),nchar(seq2))){
    if(substr(seq1,1L,i) != substr(seq2,1L,i)){
      i <- i - 1L
      break()
    }
  }
  return(substr(seq1,1L,i))
}

## Determines how many tandem repeats of seq1 are present in seq2 starting at the left
left_shared_tandem_subsequence <- function(seq1,seq2){
  stop.pos <- seq.int(from= nchar(seq1),to = (nchar(seq2) %/% nchar(seq1))*nchar(seq1),by = nchar(seq1))
  for(i in seq_along(stop.pos)){
    if(paste0(rep(seq1,i),collapse = '') != substr(seq2,1L,stop.pos[i])){
      i <- i - 1L
      break()
    }
  }
  return(substr(seq2,1L,c(0L,stop.pos)[i+1L]))
}

## Determines the smallest subequence in string by considering all subsequences that completely divide the sequence
smallest_repetitive_subsequence <- function(seq){
  len <- as.integer(nchar(seq))
  factors <- seq_len(len)[ len %% seq_len(len) == 0L]
  matches <- unlist(lapply(factors,
         function(i){
           paste0(rep(substr(seq,1,i),len/i),collapse='') == seq
         }))
  result <- substr(seq,1,min(factors[matches]))
  return(result)
}


