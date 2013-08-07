# A function for dating identified blocks of recombination.
Estimate.Ages <- function(blocks, sequence, mutation.rate=10e-8){
  stopifnot("HybRIDSdna" %in% class(sequence), is.numeric(mutation.rate))
  if("HybRIDSblock" %in% class(blocks)){
    Dates <- estimate.ages(blocks, sequence, mutation.rate)
    OutBlocks <- mergeBandD(blocks, Dates)
    OutBlocks <- list(OutBlocks, ContigNames = blocks$ContigNames)
    return(as.HybRIDSdatedBlocks(OutBlocks))
  } else {
    if("HybRIDSblockSET" %in% class(blocks)){
      cat("Now dating blocks...\n")
      DatesSet <- lapply(blocks, function(x) estimate.ages(x, sequence, mutation.rate))
      OutBlocks <- lapply(1:length(blocks), function(i) list(mergeBandD(blocks[[i]], DatesSet[[i]]), ContigNames = blocks[[i]]$ContigNames))
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




estimate.ages <- function(block, dna, mut.rate) {
  stopifnot("HybRIDSblock" %in% class(block))
  cat("Now dating blocks...\n")
  Blocks <- block[[2]]
  AllDates <- lapply(1:3, function(i) date.for.all(Blocks[[i]], dna, mut.rate, i))
  names(AllDates) <- names(block[[2]])
  return(AllDates)
}




date.for.all <- function(blocksset, Sequence, mrate, comps) {
  DatesForAllThresholds <- lapply(blocksset, function(x) date.blocks(x, Sequence, mrate, comps)) # A lapply command here because blocksset may have more than one element, 
  # depending on whether more than one suitable threshold was discovered. 
}




date.blocks <- function(blocksobj, dnaobj, mut, pair) { # blocksobj should be a dataframe
  # Check there are blocks to date!
  if(!is.character(blocksobj)){ # Checking the blocksobj is not a string of characters and does in fact contain more than one row.
    blocksobj <- as.matrix(blocksobj)
    blockAges <- matrix(nrow=nrow(blocksobj),ncol=5)
    colnames(blockAges)<-c("5%","50%","95%","BlockSize","SNPs")
    # For each significant block...
    for(i in 1:nrow(blocksobj)){
      # Pick the correct sequences for the blocks...
      blockAges[i,4] <- blocksobj[i,6]
      BlockStart <- which(colnames(dnaobj$CroppedSequence) == blocksobj[i,4])
      BlockEnd <- which(colnames(dnaobj$CroppedSequence) == blocksobj[i,5])
      # Figure out what sort of pair comparrison this comparrison is, then extract the two sequences required...
      if(pair == 1){
        Seq <- dnaobj$CroppedSequence[c(1,2),c(BlockStart:BlockEnd)]
      } else {
        if(pair == 2){
          Seq <- dnaobj$CroppedSequence[c(1,3),c(BlockStart:BlockEnd)]
        } else {
          if(pair == 3){
            Seq <- dnaobj$CroppedSequence[c(2,3),c(BlockStart:BlockEnd)]
          }
        }
      }
      # Get block length
      B <- blocksobj[i,6]
      # Make sure to remove any non-polymorphic sites and calculate the number of SNP's
      cutSeq <- as.matrix(Seq[,Seq[1,] != Seq[2,]])
      # Calculate the maximum number of SNPs
      blockAges[i,5] <- maxSNPs <- ncol(cutSeq)
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