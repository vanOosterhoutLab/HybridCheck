## Internal functions for the block dating and significant values

mergeBandD <- function( block, date ) {
  output <- lapply( 1:3, function(i) combineDatesWBlocks( block[[i]], date[[i]] ) )
  names( output ) <- names( block )
  return( output )
}

combineDatesWBlocks <- function(Bs,Ds){
  outdfs <- lapply( 1:length( Bs ), function(i) ifelse( is.character( Bs[[i]] ) || nrow( Bs[[i]] ) < 1, "NO BLOCKS DETECTED OR DATED", return( comb( Bs[[i]],Ds[[i]] ) ) ) )
  names( outdfs ) <- names( Bs )
  return( outdfs )
}

comb <- function(B,D){
  B$fiveAge <- D[,1]
  B$fiftyAge <- D[,2]
  B$ninetyfiveAge <- D[,3]
  B$SNPnum <- D[,5]
  B$PValue <- D[,6] # Incorporate the p-value into output
  B <- as.data.frame(na.omit(B))
  if( nrow( B ) == 0 ) {
    B <- "NO BLOCKS DETECTED OR DATED"
  }
  return( B )
}

binomcalc <- function(p, p0, N, B){pbinom(B,N,p)-p0}

date.blocks <- function(blocksobj, dnaobj, mut, pair, pthresh, bonfcorrect) {
  # Check there are blocks to date!
  if(!is.character(blocksobj)){ # Checking the blocksobj is not a string of characters.
    blocksobj <- as.matrix(blocksobj)
    if( bonfcorrect == TRUE ){ # Bonferoni correction of the pvalue threshold.
      pthresh <- pthresh/nrow(blocksobj)
    }
    blockAges <- matrix(nrow=nrow(blocksobj),ncol=6) #ncol was 5 before p-value accomodation.
    colnames(blockAges) <- c("5%","50%","95%","BlockSize","SNPs","p-value")
    # For each significant block...
    for(i in 1:nrow(blocksobj)){
      # Pick the correct sequences for the blocks...
      blockAges[i,4] <- blocksobj[i,6]
      BlockStart <- which(dnaobj$InformativeBp == blocksobj[i,4])
      BlockEnd <- which(dnaobj$InformativeBp == blocksobj[i,5])
      # Figure out what sort of pair comparrison this comparrison is, then extract the two sequences required...
      if(pair == 1){
        Seq <- dnaobj$InformativeSequence[c(1,2),c(BlockStart:BlockEnd)]
        wholeSequenceDist <- dist.dna(as.DNAbin(dnaobj$FullSequence[c(1,2),]), model="N")[1] # Get the raw sequence distance.
        wholeSequenceDist <- wholeSequenceDist/ncol(dnaobj$FullSequence) # Divide by the sequence length.
      } else {
        if(pair == 2){
          Seq <- dnaobj$InformativeSequence[c(1,3),c(BlockStart:BlockEnd)]
          wholeSequenceDist <- dist.dna(as.DNAbin(dnaobj$FullSequence[c(1,3),]), model="N")[1]
          wholeSequenceDist <- wholeSequenceDist/ncol(dnaobj$FullSequence)
        } else {
          if(pair == 3){
            Seq <- dnaobj$InformativeSequence[c(2,3),c(BlockStart:BlockEnd)]
            wholeSequenceDist <- dist.dna(as.DNAbin(dnaobj$FullSequence[c(2,3),]), model="N")[1]
            wholeSequenceDist <- wholeSequenceDist/ncol(dnaobj$FullSequence)
          }
        }
      }
      # Get block length
      N <- blocksobj[i,6]
      # Make sure to remove any non-polymorphic sites and calculate the number of SNP's
      cutSeq <- as.matrix(Seq[,Seq[1,] != Seq[2,]])
      # Calculate the maximum number of SNPs
      maxSNPs <- ncol( cutSeq )
      # p-value calculation is a binomial distribution, taking into account number of SNP's in the block, 
      pValue <- pbinom( maxSNPs, N, wholeSequenceDist )
      if(pValue > pthresh) next # If the block does not meet the p-value threshold, drop it and proceed to next loop iteration.
      blockAges[i,5] <- maxSNPs  
      blockAges[i,6] <- pValue
      soln5 <- uniroot(binomcalc, c(0,1), p0=0.05, B = maxSNPs, N = N)
      soln50 <- uniroot(binomcalc, c(0,1), p0=0.5, B = maxSNPs, N = N)
      soln95 <- uniroot(binomcalc, c(0,1), p0=0.95, B = maxSNPs, N = N)
      blockAges[i,3] <- soln5[["root"]]/mut
      blockAges[i,2] <- soln50[["root"]]/mut
      blockAges[i,1] <- soln95[["root"]]/mut
    }
  } else {
    blockAges <- "NO BLOCKS TO DATE"
  }
  return(blockAges)
}
