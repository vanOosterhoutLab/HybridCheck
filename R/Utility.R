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








# Define the general function - Simplify

#' @title SimplifyBlocks 
#' @description A function to simplify HybRIDSblock, HybRIDSdatedBlock, HybRIDSblockSET, and HybRIDSdatedBlocksSET objects as tables which can be exported.
#' @param x an object of class HybRIDSblock, HybRIDSdatedBlock, HybRIDSblockSET, or HybRIDSdatedBlocksSET. 
#' 
#' 
#' @rdname SimplifyBlocks
#' @export SimplifyBlocks
SimplifyBlocks <- function(x,...) UseMethod("SimplifyBlocks", x)

#' 
#' @rdname SimplifyBlocks
#' @S3method SimplifyBlocks HybRIDSblock
SimplifyBlocks.HybRIDSblock <- function(x){
  combinedframes <- data.frame(matrix(ncol=6))
  names(combinedframes) <- c("Length","Last","First","FirstBP","LastBP","ApproxBpLength")
  sequencepair <- c()
  thresholds <- c()
  for(i in 1:3){
    for(t in 1:length(x[[1]][[i]])){
      if("data.frame" %in% class(x[[2]][[i]][[t]])){
        combinedframes <- rbind(combinedframes,x[[2]][[i]][[t]])
        sequencepair <- c(sequencepair, rep(names(x[[1]])[[i]], nrow(x[[2]][[i]][[t]])))
        thresholds <- c(thresholds, rep(names(x[[2]][[i]])[[t]], nrow(x[[2]][[i]][[t]]))) # BlockInput[[1]][[i]][[t]]
      }
    }
  }
  combinedframes <- combinedframes[-1,]
  combinedframes$SequencePair <- sequencepair
  combinedframes$Threshold <- thresholds
  return(combinedframes)
}

#' 
#' @rdname SimplifyBlocks
#' @S3method SimplifyBlocks HybRIDSdatedBlocks
SimplifyBlocks.HybRIDSdatedBlocks <- function(x){
  combinedframes <- data.frame(matrix(ncol=10)) 
  names(combinedframes) <- c("Length","Last","First","FirstBP","LastBP","ApproxBpLength","fiveAge","fiftyAge","ninetyfiveAge","SNPnum")
  sequencepair <- c()
  thresholds <- c()
  for(i in 1:3){
    for(t in 1:length(x[[1]][[i]])){
      if("data.frame" %in% class(x[[1]][[i]][[t]])){
        combinedframes <- rbind(combinedframes,x[[1]][[i]][[t]])
        sequencepair <- c(sequencepair, rep(names(x[[1]])[[i]], nrow(x[[1]][[i]][[t]])))
        thresholds <- c(thresholds, rep(names(x[[1]][[i]])[[t]], nrow(x[[1]][[i]][[t]])))
      }
    }
  }
  combinedframes <- combinedframes[-1,]
  combinedframes$SequencePair <- sequencepair
  combinedframes$Threshold <- thresholds
  return(combinedframes)
}

#' 
#' @rdname SimplifyBlocks
#' @S3method SimplifyBlocks HybRIDSblockSET
SimplifyBlocks.HybRIDSblockSET  <- function(x, seq1, seq2="Individual", seq3="Pair", returnEach=T, selectionOnly=T){
  Index <- TripletIndexes(x, seq1, seq2, seq3)
  combinedframes <- data.frame(matrix(ncol=8))
  names(combinedframes) <- c("Length","Last","First","FirstBP","LastBP","ApproxBpLength","SequencePair","Threshold")
  blocksframes <- x[Index]
  sortedblocks <- lapply(blocksframes, function(x) SimplifyBlocks(x))
  for(i in 1:length(blocksframes)){
    combinedframes <- rbind(combinedframes, sortedblocks[[i]])
  }
  combinedframes <- combinedframes[-1,]
  if(returnEach==T){
    return(as.HybRIDSdatedBlocksTableSET(sortedblocks))
  } else {
    # We need to make sure that there are no duplicated in the DF if there are we can try to count them.
    if(selectionOnly == T && seq3=="Pair"){
      combinedframes <- combinedframes[combinedframes$SequencePair == paste(seq1, seq2, sep=":"),]
    } else {
      if(selectionOnly == T && seq3=="Pair" && seq2=="Individual"){
        combinedframes <- combinedframes[unlist(lapply(strsplit(combinedframes$SequencePair, ":"), function(x) seq1 %in% x)),]
      }
    }
    # Make sure there are no duplicates
    dups <- duplicated(combinedframes)
    if(!any(dups)) return(as.HybRIDSblockTable(combinedframes))
  }
  return(as.HybRIDSblockTable(combinedframes))
  # This bit of the code may seem silly but I have plans for it later.
}

#' 
#' @rdname SimplifyBlocks
#' @S3method SimplifyBlocks HybRIDSdatedBlocksSET
SimplifyBlocks.HybRIDSdatedBlocksSET  <- function(x, seq1, seq2="Individual", seq3="Pair", returnEach=T, selectionOnly=T){
  Index <- TripletIndexes(x, seq1, seq2, seq3)
  combinedframes <- data.frame(matrix(ncol=12))
  names(combinedframes) <- c("Length","Last","First","FirstBP","LastBP","ApproxBpLength","fiveAge","fiftyAge","ninetyfiveAge","SNPnum","SequencePair","Threshold")
  blocksframes <- x[Index]
  sortedblocks <- lapply(blocksframes, function(x) SimplifyBlocks(x))
  for(i in 1:length(blocksframes)){
    combinedframes <- rbind(combinedframes, sortedblocks[[i]])
  }
  combinedframes <- combinedframes[-1,]
  if(returnEach==T){
    return(as.HybRIDSblockTableSET(sortedblocks))
  } else {
    # We need to make sure that there are no duplicated in the DF if there are we can try to count them.
    if(selectionOnly == T && seq3=="Pair"){
      combinedframes <- combinedframes[combinedframes$SequencePair == paste(seq1, seq2, sep=":"),]
    } else {
      if(selectionOnly == T && seq3=="Pair" && seq2=="Individual"){
        combinedframes <- combinedframes[unlist(lapply(strsplit(combinedframes$SequencePair, ":"), function(x) seq1 %in% x)),]
      }
    }
    # Make sure there are no duplicates
    dups <- duplicated(combinedframes)
    if(!any(dups)) return(as.HybRIDSblockTable(combinedframes))
  }
  return(as.HybRIDSblockTable(combinedframes))
  # This bit of the code may seem silly but I have plans for it later.
}