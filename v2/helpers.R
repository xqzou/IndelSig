################################################################################
## Helper functions for seq processing
################################################################################
## Determines the shared subsequence between two sequences starting at the left
left_shared_subsequence <- function(seq1,seq2){
 comparison <- unlist(lapply(1:min(nchar(seq1),nchar(seq2)),
                            function(i){
                              substr(seq1,1,i) == substr(seq2,1,i)
                              }))
 comparison <- substr(seq1,1,sum(comparison))
 return(comparison)
}

## Determines how many tandem repeats of seq1 are present in seq2 starting at the left
left_shared_tandem_subsequence <- function(seq1,seq2){
  stop.pos <- seq(from= nchar(seq1),to = (nchar(seq2) %/% nchar(seq1))*nchar(seq1),by = nchar(seq1))
  comparison <- unlist(lapply(stop.pos,
         function(i){
           paste0(rep(seq1,i %/% nchar(seq1)),collapse='')==substr(seq2,1,i)
         }))
  comparison <- substr(seq2,1,c(0,stop.pos)[sum(comparison)+1])
  return(comparison)
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


