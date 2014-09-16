# Utilities and Functions for ease of programmatically interacting with HybRIDS objects.
# (And to make coding the GUI more concise!)

#' Subset a HybRIDSdna object and create a new HybRIDSdna object. 
#' 
#' Create a new HybRIDSdna object from a selection of sequences present in another HybRIDSdna object.
#' 
#' This is a convienience utility function to save the user from writing a chunk of code playing the names and indexing the dataset.
#' It is also an important function for the sequence selector in the HybRIDS GUI.
#' @param x A HybRIDSdna object
#' @param sequences Either a vector of sequence numbers, or sequence names, the function automaticaly detects which.
#' @export
SubSeq <- function(x, sequences){
  if(!"HybRIDSdna" %in% class(x)) stop("You need to enter a valid HybRIDSdna object for argument x!")
  if(!"character" %in% class(sequences) && !"numeric" %in% class(sequences)) stop("Argument: sequences needs to be a vector of numbers or character strings for either sequence number or sequence names")
  if("character" %in% class(sequences)) {
    seqnumbers <- which(x$ContigNames %in% svalue(sequences))
  } else {
    seqnumbers <- sequences
  }
  output <- x
  output$Sequence <- output$Sequence[seqnumbers,]
  output$CroppedSequence <- output$CroppedSequence[seqnumbers,]
  output$ContigNames <- output$ContigNames[seqnumbers]
  combs <- combn(1:nrow(output$Sequence),3, simplify=FALSE)
  output$Combinations <- combs
  return(output)
}

# Define Generic Method for Subsetting Triplets based on sequence names...

#' TripletIndexes
#' 
#' A function to get the numeric index of the triplet containing the sequence names.
#' 
#' @param SSSet An object that is either a HybRIDSseqsimSET, HybRIDSblockSET, or HybRIDSdatedBlocksSET
#' @param sequence1 The name of the first sequence.
#' @param sequence2 The name os the second sequence.
#' @param sequence3 The same of the third sequence.
#' @details If a sequence name is provided for parameters sequence1, sequence2, and sequence3, then a single index corresponding to the triplet in which all three sequences occur is returned.
#' Otherwise if sequence1 and 2 are provided but sequence3 is not (so is left as "Pair"), several indicies are returned which are triplets containing the pair specified. 
#' Furthuremore if sequence1 is provided (so sequence3 left as "Pair" and sequence 2 left as "Individual"),  several indicies are returned which are triplets containing the single seqence specified
TripletIndexes <- function(SSSet, sequence1, sequence2="Individual", sequence3="Pair"){
  if(length(unique(c(sequence1, sequence2, sequence3))) == 1){
    stop("You can't set all the sequence names to the same!")
  }
  if(sequence3 == "Pair" && sequence2 == "Individual"){
      outindex <- which(unlist(lapply(SSSet, function(x) sequence1 %in% x$ContigNames)))
  } else {
    if(sequence3 == "Pair" && sequence2 != "Individual"){
      outindex <- which(unlist(lapply(SSSet, function(x) all(c(sequence1, sequence2) %in% x$ContigNames))))
    } else {
      if(sequence3 != "Pair" && sequence2 != "Individual"){
        outindex <- which(unlist(lapply(SSSet, function(x) all(c(sequence1, sequence2, sequence3) %in% x$ContigNames))))
      } else {
        if(sequence3 != "Pair" && sequence2 == "Individual"){
          stop("Can't do this - You've specified 'Individual' but provided two sequence names...")
        }
      }
    }
  }
  return(outindex)
}