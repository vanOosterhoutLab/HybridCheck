#' Function that estimates the ages of putative recombination blocks
#' 
#' This function take a HybRIDSblock or HybRIDSblockSET and will estimate the ages of the blocks detected.
#' 
#' @param blocks An object of class HybRIDSblock or HybRIDSblockSET
#' @param sequence The object of class HybRIDSdna that was used to generate the object given as the 'blocks' parameter.
#' @export
Estimate.Ages <- function(blocks, sequence, mutation.rate = 10e-8, requiredP = 0.005){
  stopifnot("HybRIDSdna" %in% class(sequence), is.numeric(mutation.rate))
  if("HybRIDSblock" %in% class(blocks)){
    Dates <- estimate.ages(blocks, sequence, mutation.rate, requiredP)
    OutBlocks <- mergeBandD(blocks, Dates)
    OutBlocks <- list(OutBlocks, ContigNames = blocks$ContigNames)
    return(as.HybRIDSdatedBlocks(OutBlocks))
  } else {
    if("HybRIDSblockSET" %in% class(blocks)){
      cat("Now dating blocks...\n")
      DatesSet <- lapply(blocks, function(x) estimate.ages(x, sequence, mutation.rate, requiredP))
      OutBlocks <- lapply(1:length(blocks), function(i) as.HybRIDSdatedBlocks(list(mergeBandD(blocks[[i]], DatesSet[[i]]), ContigNames = blocks[[i]]$ContigNames))) # Put as.HybRIDSdatedBlocks here to make sure each element of the set is correct type...
      return(as.HybRIDSdatedBlocksSET(OutBlocks))
    } else {
      stop("Blocks input was not of type HybRIDSblock or HybRIDSblockSET")
    }
  }
}


mergeBandD <- function(block, date){
  combineDatesWBlocks <- function(Bs,Ds){
    comb <- function(B,D){
      B$fiveAge <- D[,1]
      B$fiftyAge <- D[,2]
      B$ninetyfiveAge <- D[,3]
      B$SNPnum <- D[,5]
      B$PValue <- D[,6] # Incorporate the p-value into output
      B <- as.data.frame(na.omit(B))
      if(nrow(B) == 0){
        B <- "NO BLOCKS DETECTED OR DATED"
      }
      return(B)
    }
    outdfs <- lapply(1:length(Bs), function(i) ifelse(is.character(Bs[[i]]) || nrow(Bs[[i]]) < 1, "NO BLOCKS DETECTED OR DATED", return(comb(Bs[[i]],Ds[[i]]))))
    names(outdfs) <- names(Bs)
    return(outdfs)
  }
  output <- lapply(1:3, function(i) combineDatesWBlocks(block[[2]][[i]], date[[i]]))
  names(output) <- names(block[[2]])
  return(output)
}




estimate.ages <- function(block, dna, mut.rate, rp) {
  stopifnot("HybRIDSblock" %in% class(block))
  cat("Now dating blocks...\n")
  Blocks <- block[[2]]
  AllDates <- lapply(1:3, function(i) date.for.all(Blocks[[i]], dna, mut.rate, i, rp))
  names(AllDates) <- names(block[[2]])
  return(AllDates)
}




date.for.all <- function(blocksset, Sequence, mrate, comps, reqp) {
  DatesForAllThresholds <- lapply(blocksset, function(x) date.blocks(x, Sequence, mrate, comps, reqp)) # A lapply command here because blocksset may have more than one element, 
  # depending on whether more than one suitable threshold was discovered. 
}




date.blocks <- function(blocksobj, dnaobj, mut, pair, pthresh) {
  # Check there are blocks to date!
  if(!is.character(blocksobj)){ # Checking the blocksobj is not a string of characters and does in fact contain more than one row.
    blocksobj <- as.matrix(blocksobj)
    blockAges <- matrix(nrow=nrow(blocksobj),ncol=6) #ncol was 5 before p-value accomodation.
    colnames(blockAges) <- c("5%","50%","95%","BlockSize","SNPs","p-value")
    # For each significant block...
    for(i in 1:nrow(blocksobj)){
      # Pick the correct sequences for the blocks...
      blockAges[i,4] <- blocksobj[i,6]
      BlockStart <- which(colnames(dnaobj$CroppedSequence) == blocksobj[i,4])
      BlockEnd <- which(colnames(dnaobj$CroppedSequence) == blocksobj[i,5])
      # Figure out what sort of pair comparrison this comparrison is, then extract the two sequences required...
      if(pair == 1){
        Seq <- dnaobj$CroppedSequence[c(1,2),c(BlockStart:BlockEnd)]
        wholeSequenceDist <- dist.dna(as.DNAbin(dnaobj$Sequence[c(1,2),]), model="N")[1] # Get the raw sequence distance.
        wholeSequenceDist <- wholeSequenceDist/ncol(dnaobj$Sequence) # Divide by the sequence length.
        #wholeSequenceDist <- wholeSequenceDist*100 # Multiply by 100 to get divergence of two sequences as a percentage.
      } else {
        if(pair == 2){
          Seq <- dnaobj$CroppedSequence[c(1,3),c(BlockStart:BlockEnd)]
          wholeSequenceDist <- dist.dna(as.DNAbin(dnaobj$Sequence[c(1,3),]), model="N")[1]
          wholeSequenceDist <- wholeSequenceDist/ncol(dnaobj$Sequence)
          #wholeSequenceDist <- wholeSequenceDist*100
        } else {
          if(pair == 3){
            Seq <- dnaobj$CroppedSequence[c(2,3),c(BlockStart:BlockEnd)]
            wholeSequenceDist <- dist.dna(as.DNAbin(dnaobj$Sequence[c(2,3),]), model="N")[1]
            wholeSequenceDist <- wholeSequenceDist/ncol(dnaobj$Sequence)
            #wholeSequenceDist <- wholeSequenceDist*100
          }
        }
      }
      # Get block length
      B <- blocksobj[i,6]
      # Make sure to remove any non-polymorphic sites and calculate the number of SNP's
      cutSeq <- as.matrix(Seq[,Seq[1,] != Seq[2,]])
      # Calculate the maximum number of SNPs
      maxSNPs <- ncol(cutSeq)
      
      # p-value calculation is a binomial distribution, taking into account number of SNP's in the block, 
      pValue <- pbinom(maxSNPs, B, wholeSequenceDist)
      if(pValue > pthresh) next # If the block does not meet the p-value threshold, drop it and proceed to next loop iteration. 
      
      blockAges[i,5] <- maxSNPs  
      blockAges[i,6] <- pValue
      t <- 0
      continue <- T
      flag5 <- F
      flag50 <- F
      flag95 <- F
      while(continue==T){
        p <- mut*t
        sumprobs <- pbinom(maxSNPs,B,p)*1000
        if(sumprobs<=50L && flag5==F){
          blockAges[i,3] <- t
          flag5 <- T
        } else {
          if(sumprobs<=500L && flag50==F){
            blockAges[i,2] <- t
            flag50 <- T
          } else {
            if(sumprobs<=950L && flag95==F){
              blockAges[i,1] <- t
              flag95 <- T
            }
          }
        }
        if(all(c(flag5,flag50,flag95) == TRUE)){
          continue <- F
        }
        t <- t+100
      }
    }
  } else {
    blockAges <- "NO BLOCKS TO DATE"
  }
  return(blockAges)
}