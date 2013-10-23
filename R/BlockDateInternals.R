## Internal functions for the block dating and significant values.

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

date.blocks <- function(blocksobj, dnaobj, mut, pair, pthresh, bonfcorrect, danyway) {
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
      #Extract the two sequences required...
      Seq <- dnaobj$InformativeSequence[pair,c(BlockStart:BlockEnd)]
      wholeSequenceDist <- (dist.dna(as.DNAbin(dnaobj$FullSequence[pair,]), model="N")[1])/ncol(dnaobj$FullSequence) # Get the raw sequence distance.
      # Get block length
      N <- blocksobj[i,6]
      # Make sure to remove any non-polymorphic sites and calculate the number of SNP's
      maxSNPs <- length(seg.sites(as.DNAbin(Seq)))
      # Probability of observing this many SNPs in the block given the divergence of the whole sequence divergence of the sequences. 
      pValue <- pbinom( maxSNPs, N, wholeSequenceDist )
      if( pValue > pthresh && danyway == FALSE ) {
        next # If the block does not meet the p-value threshold, drop it and proceed to next loop iteration.
      } else {
        blockAges[i,5] <- maxSNPs  
        blockAges[i,6] <- pValue
        soln5 <- uniroot(binomcalc, c(0,1), p0=0.05, B = maxSNPs, N = N)
        soln50 <- uniroot(binomcalc, c(0,1), p0=0.5, B = maxSNPs, N = N)
        soln95 <- uniroot(binomcalc, c(0,1), p0=0.95, B = maxSNPs, N = N)
        blockAges[i,3] <- soln5[["root"]]/mut
        blockAges[i,2] <- soln50[["root"]]/mut
        blockAges[i,1] <- soln95[["root"]]/mut
      }
    }
  } else {
    blockAges <- "NO BLOCKS TO DATE"
  }
  return(blockAges)
}
