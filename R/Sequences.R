### Manipulation of DNA sequences.

## General functions that may be useful outside of their internal use in
## HybridCheck.

# Return alignment of only the selected base positions.
extractBases <- function(dna, indexes){
  if(!is(dna, "DNAStringSet")){stop("Argument dna needs to be of class DNAStringSet")}
  if(!is.numeric(indexes)){stop("Argument indexes needs to be of class Numeric")}
  i <- rep.int(IntegerList(indexes), length(dna))
  newSeq <- dna[i]
  return(newSeq)
}


# Find inforative sites in DNAStringSet and return a DNAStringSet
# containing only the informative sites.
cutToInformativeSites <- function(sequences, sequenceInfo){
  conMat <- consensusMatrix(sequences)
  sequenceInfo$InformativeUsed <- 
    which(
      (colSums(conMat[c(5:15, 17, 18),]) == 0) &
        (colSums(conMat != 0) > 1)
    )
  sequenceInfo$InformativeUsedLength <- 
    length(sequenceInfo$InformativeUsed)
  return(extractBases(sequences, sequenceInfo$InformativeUsed))
}