## Internal functions for the block dating and significant values.

mergeBandD <- function(block, date) {
  output <- lapply(1:3, function(i) combineDatesWBlocks(block[[i]], date[[i]]))
  names(output) <- names(block)
  return(output)
}

combineDatesWBlocks <- function(Bs,Ds){
  outdfs <- lapply(1:length(Bs), function(i) ifelse(is.character(Bs[[i]]) || nrow(Bs[[i]]) < 1, "NO BLOCKS DETECTED OR DATED", return(comb(Bs[[i]], Ds[[i]]))))
  names(outdfs) <- names(Bs)
  return(outdfs)
}

comb <- function(B,D){
  B$fiveAge <- D[,1]
  B$fiftyAge <- D[,2]
  B$ninetyfiveAge <- D[,3]
  B$SNPnum <- D[,5]
  B$PValue <- D[,6] # Incorporate the p-value into output
  B$PThresh <- D[,7]
  B <- as.data.frame(na.omit(B))
  if(nrow(B) == 0) {
    B <- "NO BLOCKS DETECTED OR DATED"
  }
  return(B)
}

# #########################################################################################################################
# 
mutation.correction <- function(ObservedMutations, SequenceSize){
  ObservedRatio <- ObservedMutations/SequenceSize
  ActualRatio <- seq(from=0, to=2, by=0.0001)
  OutputRatio <- (0.992582633 * ActualRatio) - (0.605566567 * (ActualRatio^2) ) + (0.166571989 * (ActualRatio^3) )
  ResultRatio <- ActualRatio[which(OutputRatio >= ObservedRatio)[1]]
  #Output <- ceiling(ResultRatio * SequenceSize)
  return(ResultRatio)
}
# 
# ActualRatio <- seq(from=0, to=1, by=0.0001)
# OutputRatio <- (0.9980 * ActualRatio) - (0.6170 * (ActualRatio^2) ) + (0.1704 * (ActualRatio^3) )
# ResultRatio <- ActualRatio[which(OutputRatio >= ObservedRatio)[1]]
# 
# OutputRatio <- (0.992582633 * ActualRatio) - (0.605566567 * (ActualRatio^2) ) + (0.166571989 * (ActualRatio^3) )
# 
# 
# results_10000_03$testresult[which(is.na(results_10000_03$testresult))] <- apply(results_1000_04[which(is.na(results_10000_03$testresult)),], 1, function(x) mutation.correction(x[1], x[4]))
# 
# lm(results_1000_03$ObservedSNPsRatio ~ results_1000_03$ActualSNPsRatio + I(results_1000_03$ActualSNPsRatio^2) + I(results_1000_03$ActualSNPsRatio^3))
# 
# ActualRatio <- seq(from=0, to=2, by=0.0001)
# OutputRatio <- (0.992582633 * ActualRatio) - (0.605566567 * (ActualRatio^2) ) + (0.166571989 * (ActualRatio^3) )
# 
# model <- lm(results_10000_03$ActualSNPsRatio ~ results_1000_03$ObservedSNPsRatio + I(results_1000_03$ObservedSNPsRatio^2) + I(results_1000_03$ObservedSNPsRatio^3))
# x <- seq(from=0, to=0.7, by=0.1)
# plot(results_10000_03$ObservedSNPsRatio, results_10000_03$ActualSNPsRatio)
# points(predict(model), col="red")
# 
# results_5000_02$testresult <- apply(results_5000_02, 1, function(x) mutation.correction(x[1], x[4]))
# results_100_02$testresult <- apply(results_100_02, 1, function(x) mutation.correction(x[1], x[4]))
# results_1000_03$testresult <- apply(results_1000_03, 1, function(x) mutation.correction(x[1], x[4]))
# results_10000_03$testresult <- apply(results_10000_03, 1, function(x) mutation.correction(x[1], x[4]))
# results_100_03$testresult <- apply(results_100_03, 1, function(x) mutation.correction(x[1], x[4]))
# results_5000_03$testresult <- apply(results_5000_03, 1, function(x) mutation.correction(x[1], x[4]))
# results_10000_04$testresult <- apply(results_10000_04, 1, function(x) mutation.correction(x[1], x[4]))
# results_1000_04$testresult <- apply(results_1000_04, 1, function(x) mutation.correction(x[1], x[4]))  
# results_5000_04$testresult <- apply(results_5000_04, 1, function(x) mutation.correction(x[1], x[4]))
# 
# save(results_1000_03,file="results_1000_03_highreps")
# save(results_1000_03,file="results_1000_03_highreps")
# save(results_10000_03,file="results_10000_03_highreps")
# save(results_1000_04,file="results_1000_04_highreps")
# save(results_10000_04,file="results_10000_04_highreps")
# save(results_10000_05,file="results_10000_05_highreps")
# save(results_100_03,file="results_100_03_highreps")
# save(results_100_02,file="results_100_02_highreps")
# save(results_5000_04,file="results_5000_04_highreps")
# save(results_5000_03,file="results_5000_03_highreps")
# save(results_5000_02,file="results_5000_02_highreps")
# 
# ##########################################################################################################


binomcalc <- function(p, p0, N, B){pbinom(B,N,p)-p0}

date.blocks <- function(blocksobj, dnaobj, mut, pair, pthresh, bonfcorrect, danyway) {
  # Check there are blocks to date!
  if(!is.character(blocksobj)){ # Checking the blocksobj is not a string of characters.
    blocksobj <- as.matrix(blocksobj)
    if( bonfcorrect == TRUE ){ # Bonferoni correction of the pvalue threshold.
      pthresh <- pthresh/nrow(blocksobj)
    }
    blockAges <- matrix(nrow=nrow(blocksobj),ncol=7) #ncol was 6 before p-thresh accomodation.
    colnames(blockAges) <- c("5%","50%","95%","BlockSize","SNPs","p-value", "p-threshold")
    # For each significant block...
    for(i in 1:nrow(blocksobj)){
      # Pick the correct sequences for the blocks...
      blockAges[i,4] <- blocksobj[i,"ApproxBpLength"]
      BlockStart <- which(dnaobj$FullBp == blocksobj[i,"FirstBP"]) # Was Informativebp
      BlockEnd <- which(dnaobj$FullBp == blocksobj[i,"LastBP"]) # Was informative bp
      #Extract the two sequences required...
      Seq <- dnaobj$FullSequence[pair,c(BlockStart:BlockEnd)]
      wholeSequenceDist <- (dist.dna(as.DNAbin(dnaobj$FullSequence[pair,]), model="raw")[1])
      # Get block length
      N <- ncol(Seq)
      # Make sure to remove any non-polymorphic sites and calculate the number of SNP's
      maxSNPs <- dist.dna(as.DNAbin(Seq), model="N")[1]
      # Probability of observing this many SNPs in the block given the divergence of the whole sequence divergence of the sequences. 
      pValue <- pbinom( maxSNPs, N, wholeSequenceDist )
      #maxSNPs <- mutation.correction(maxSNPs, N)
      if( pValue > pthresh && danyway == FALSE ) {
        next # If the block does not meet the p-value threshold, drop it and proceed to next loop iteration.
      } else {
        if(N == maxSNPs || pValue == 1) next
        blockAges[i,5] <- maxSNPs  
        blockAges[i,6] <- pValue
        soln5 <- uniroot(binomcalc, c(0,1), p0=0.05, B=maxSNPs, N=N)
        soln50 <- uniroot(binomcalc, c(0,1), p0=0.5, B=maxSNPs, N=N)
        soln95 <- uniroot(binomcalc, c(0,1), p0=0.95, B=maxSNPs, N=N)
        blockAges[i,3] <- round(soln5[["root"]]/(2*mut))
        blockAges[i,2] <- round(soln50[["root"]]/(2*mut))
        blockAges[i,1] <- round(soln95[["root"]]/(2*mut))
        blockAges[,7] <- pthresh
      }
    }
  } else {
    blockAges <- "NO BLOCKS TO DATE"
  }
  return(blockAges)
}
