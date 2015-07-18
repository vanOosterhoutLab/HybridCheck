# Resolving the ambigious base positions in DNA sequences.
# Ben J. Ward, 2015. 


# Functions for HCseq.

## Check population selection are names or integers.
popIntegersToNames <- function(popSel, seqNames){
  pops <- lapply(popSel, function(x){
    if(is.integer(x)){
      return(seqNames[x])
    } else {
      if(is.character(x)){
        return(x)
      } else {
        stop("Need to provide a list of groups of sequence names or integers representing sequence numbers.")
      }
    }
  })
  if(any(table(unlist(pops)) > 1)){
    stop("Entered a sequence name or number in more than one group.")
  }
  if(any(!unlist(lapply(pops, function(x) all(x %in% seqNames))))){
    stop("Some sequences specified in the populations are not in the sequence data.")
  }
  return(pops)
}


compDist <- function(popPairs, seqInPop, distMat){
  distances <- lapply(popPairs, function(y){
    grid <- expand.grid(seqInPop[y])
    ds <- apply(grid, 1, function(z){distMat[z[1], z[2]]})
    return(sum(ds) / nrow(grid))
  })
  out <- data.frame(do.call(rbind, popPairs))
  colnames(out) <- c("OTU1", "OTU2")
  out$OTU1 <- as.character(out$OTU1)
  out$OTU2 <- as.character(out$OTU2)
  out$dist <- unlist(distances)
  return(out)
}

calculateStats <- function(counts.all, biSites.all, slice1, slice2, slice3, slice4){
  # Allocate space for Cabba and Cbaba.
  ABBA <- vector(mode = "numeric", length = length(biSites.all))
  BABA <- vector(mode = "numeric", length = length(biSites.all))
  # Allocate space for maxABBA_D.
  maxABBA_23 <- vector(mode = "numeric", length = length(biSites.all))
  maxBABA_23 <- vector(mode = "numeric", length = length(biSites.all))
  maxABBA_13 <- vector(mode = "numeric", length = length(biSites.all))
  maxBABA_13 <- vector(mode = "numeric", length = length(biSites.all))
  for(i in 1:length(biSites.all)){
    i.bi.counts <- counts.all[, biSites.all[i]]
    alleles <- names(which(i.bi.counts > 0))
    anc <- names(which(slice4$countsBi[, i] != 0))
    if(length(anc) != 1){
      anc <- names(which(i.bi.counts == max(i.bi.counts)))[1]
    }
    derived <- alleles[which(alleles != anc)]
    P1df <- slice1$countsBi[derived, i] / sum(slice1$countsBi[, i])
    P2df <- slice2$countsBi[derived, i] / sum(slice2$countsBi[, i])
    P3df <- slice3$countsBi[derived, i] / sum(slice3$countsBi[, i])
    P4df <- slice4$countsBi[derived, i] / sum(slice4$countsBi[, i])
    ABBA[i] <- (1 - P1df) * P2df * P3df * (1 - P4df)
    BABA[i] <- P1df * (1 - P2df) * P3df * (1 - P4df)
    if(!is.na(P3df) & !is.na(P2df) & P3df >= P2df){
      maxABBA_23[i] <- (1 - P1df) * P3df * P3df * (1 - P4df)
      maxBABA_23[i] <- P1df * (1 - P3df) * P3df * (1 - P4df)
    } else {
      maxABBA_23[i] <- (1 - P1df) * P2df * P2df * (1 - P4df)
      maxBABA_23[i] <- P1df * (1 - P2df) * P2df * (1 - P4df)
    }
    if(!is.na(P3df) & !is.na(P1df) & P3df >= P1df){
      maxABBA_13[i] <- (1 - P3df) * P2df * P3df * (1 - P4df)
      maxBABA_13[i] <- P3df * (1 - P2df) * P3df * (1 - P4df)
    } else {
      maxABBA_13[i] <- (1 - P1df) * P2df * P1df * (1 - P4df)
      maxBABA_13[i] <- P1df * (1 - P2df) * P1df * (1 - P4df)
    }
  }
  numABBA <- length(which(ABBA > BABA))
  numBABA <- length(which(ABBA < BABA))
  if(numABBA < numBABA){
    binomialP <- pbinom(numABBA, (numABBA + numBABA), .5, lower.tail = TRUE)
  } else {
    binomialP <- pbinom(numBABA, (numABBA + numBABA), .5, lower.tail = TRUE)
  }
  out <- data.frame(numABBA = numABBA, numBABA = numBABA, numBinomialP = binomialP,
                    ABBA = sum(ABBA), BABA = sum(BABA), maxABBA_23 = sum(maxABBA_23),
                    maxBABA_23 = sum(maxBABA_23), maxABBA_13 = sum(maxABBA_13),
                    maxBABA_13 = sum(maxBABA_13))
  return(out)
}


populationSlice <- function(popSeqs, biSites){
  counts <- consensusMatrix(popSeqs)
  alleles <- colSums(counts != 0)
  S <- sum(alleles > 1)
  return(list(counts = counts, countsBi = counts[, biSites], alleles = alleles, S = S))
}

calculateDandFd <- function(aln, pops){
  # State counts at each site.
  counts.all <- consensusMatrix(aln)
  # Calculate the number of alleles at each site in the alignment.
  message("\t\t- Calculating D and Fd for jack-knife block.")
  num.alleles.all <- colSums(counts.all != 0)
  # Find which sites are bi-allelic.
  bi.all <- which(num.alleles.all == 2)
  # Find which ones are variable.
  var.sites.all <- which(num.alleles.all > 1)
  # Find the sites containing uncertain bases.
  uncertain.all <- which(colSums(counts.all[c(5:15, 17:18), ] != 0) >= 1)
  # Make sure that bi.all and var.sites.all does not include sites with uncertain bases.
  bi.all <- bi.all[-which(bi.all %in% uncertain.all)]
  var.sites.all <- var.sites.all[-which(var.sites.all %in% uncertain.all)]
  # Make a version of the sequence alignment which only includes variable sites.
  aln.var <- subsetSequence(aln, var.sites.all)
  # Get the number of alleles at those variable sites.
  alleles.var <- num.alleles.all[var.sites.all]
  s.all <- length(alleles.var)
  # Find which of the variable sites are bi-allelic.
  bi.var <- which(alleles.var == 2)
  # Population slices - refer to function's description.
  p1Slice <- populationSlice(aln.var[pops[[1]]], bi.var)
  p2Slice <- populationSlice(aln.var[pops[[2]]], bi.var)
  p3Slice <- populationSlice(aln.var[pops[[3]]], bi.var)
  p4Slice <- populationSlice(aln.var[pops[[4]]], bi.var)
  return(data.frame(S = s.all, P1_S = p1Slice$S, P2_S = p2Slice$S,
                    P3_S = p3Slice$S, P4_S = p4Slice$S, 
                    calculateStats(counts.all, bi.all,
                                   p1Slice, p2Slice, p3Slice, p4Slice)))
}

#' @title fourTaxonTest
#' @name fourTaxonTest
#' @description Computes the Four Taxon Test for a given set of P1, P2, P3, and P4.
#' Patterson's D is calculated, as is Fd from the paper Martin et al. (2014).
#' Fd is calculated for the scenario of complete introgression between P1 and P3,
#' and also for the scenario of complete introgression between P2 and P3.
#' The jack-knife implementation and computation of Z scores is based on that used in the software package
#' ANGSD.
fourTaxonTest <- function(dna, fttRecord, numBlocks, lengthOfBlocks){
  # The jack-knife implementation and computation of Z score has been based on that
  # used in the software ANGSD.
  message(" - Running a Four Taxon Test for the four populations:\n    ",
          paste(fttRecord$P1, fttRecord$P2, 
                fttRecord$P3, fttRecord$A, collapse = ", "))
  message("\t- Figureing out the size and number of jack-knife blocks.")
  if(!is.numeric(numBlocks) && !is.numeric(lengthOfBlocks)){
    stop("Invalid input - no number of blocks or size of blocks provided.")
  }
  if(is.numeric(numBlocks) && is.numeric(lengthOfBlocks)){
    stop("Provide the number of blocks to divide sequence alignment into, OR the proposed size of the blocks, not both.")
  }
  dnaLen <- dna$getFullLength()
  if(is.numeric(lengthOfBlocks) && !is.numeric(numBlocks)){
    fttRecord$numBlocks <- as.integer(floor(dnaLen / lengthOfBlocks))
  } else {
    fttRecord$numBlocks <- numBlocks
  }
  fttRecord$blockLength <- as.integer(floor(dnaLen / fttRecord$numBlocks))
  blockStart <- seq(from = 1, to = dnaLen, by = fttRecord$blockLength)
  blockEnd <- seq(from = fttRecord$blockLength, to = dnaLen, by = fttRecord$blockLength)
  if(length(blockStart) > length(blockEnd)){
    blockStart <- blockStart[1:length(blockStart) - abs(length(blockStart) - length(blockEnd))]
  }
  # Calculation of ABBA and BABA from segments of the alignment:
  results <- data.frame(BlockStart = blockStart, BlockEnd = blockEnd)
  
  # SORT THIS OUT - IT IS INNEFICIENT TO SUBSET SEQUENCES LIKE THIS
  message("\t- Subsetting sequences for jack-knife blocks...")
  
  #blocks <- apply(results, 1, function(x){subsetSequence(dna$FullSequence, x[1]:x[2])})
  blocks <- apply(results, 1, function(x){
    subseq(dna$FullSequence, start = x[1], end = x[2])
    })
  
  message("\t- Calculating statistics needed for D and Fd, for each jack-knife block.")
  
  blocksStats <- do.call(rbind, lapply(blocks, function(x){
    calculateDandFd(x, dna$Populations[c(fttRecord$P1, fttRecord$P2, fttRecord$P3, fttRecord$A)])}))
  
  message("\t- Calculating full set of Four Taxon Test statistics for all jackknife blocks.")
  # Calculation of stats for each jackknife segment:
  # Calculation of Observed S and S for complete introgression scenarios between P1:P3, and P2:P3.
  blocksStats$S_1234 <- blocksStats$ABBA - blocksStats$BABA
  blocksStats$S_1DD4 <- blocksStats$maxABBA_23 - blocksStats$maxBABA_23
  blocksStats$S_D2D4 <- blocksStats$maxABBA_13 - blocksStats$maxBABA_13
  # Calculation of the observed pattersons D value per segment.
  blocksStats$abbaBabaSum <- blocksStats$ABBA + blocksStats$BABA
  blocksStats$D <- blocksStats$S_1234 / blocksStats$abbaBabaSum
  # Calculation of the Fd values for the two complete introgression scenarios per segment.
  blocksStats$Fd_1DD4 <- blocksStats$S_1234 / blocksStats$S_1DD4
  blocksStats$Fd_D2D4 <- blocksStats$S_1234 / blocksStats$S_D2D4
  blocksStats$Fd_1DD4_D0 <- as.numeric(blocksStats$Fd_1DD4)
  blocksStats$Fd_1DD4_D0[c(which(is.na(blocksStats$D)), which(blocksStats$D < 0))] <- NA
  blocksStats$Fd_D2D4_D0 <- as.numeric(blocksStats$Fd_D2D4)
  blocksStats$Fd_D2D4_D0[c(which(is.na(blocksStats$D)), which(blocksStats$D > 0))] <- NA
  
  # Jack-knifing based on ANGSD implementation.
  # Number of blocks.
  nBlocks <- nrow(blocksStats)
  
  # Calculates D. 
  statCalc <- function(x){sum(x[, 1])/sum(x[, 2])}
  
  # Global Estimates.
  message("\t- Calculating Global Four Taxon Test statistics.")
  globalEstimate_D <- statCalc(blocksStats[, c("S_1234", "abbaBabaSum")])
  globalEstimate_Fd_1DD4 <- statCalc(blocksStats[, c("S_1234", "S_1DD4")])
  globalEstimate_Fd_D2D4 <- statCalc(blocksStats[, c("S_1234", "S_D2D4")])
  
  # Pseudoestimates:
  
  message("\t- Making pseudoestimates as part of jack-knifing process.")
  # Allocate space for the results.
  blocksStats$pseudoD <- blocksStats$pseudoFd_1DD4 <- blocksStats$pseudoFd_D2D4 <- rep(0, nBlocks)
  
  # blockFraction shows the proportion of all ABBA and BABA sites which are located in a given jack-knife block.
  # It quantifies the influence the block has over the global result.
  blocksStats$blockFraction <- blocksStats$abbaBabaSum / sum(blocksStats$abbaBabaSum)
  
  # Loop over and count up the pseudoestimates of D and Fd:
  for(i in 1:nBlocks){
    blocksStats$pseudoD[i] <- statCalc(blocksStats[-i, c("S_1234", "abbaBabaSum")])
    blocksStats$pseudoFd_1DD4[i] <- statCalc(blocksStats[-i, c("S_1234", "S_1DD4")])
    blocksStats$pseudoFd_D2D4[i] <- statCalc(blocksStats[-i, c("S_1234", "S_D2D4")])
  }
  
  # Calculate the inverse of the block fraction - this is how much the pseudo values should be multiplied by
  # to correct for their influence given the amount of two-state sites the jack-knife block accounts for.
  blocksStats$invBlockFraction <- 1 - blocksStats$blockFraction
  
  # Calculate scaled pseudo values.
  blocksStats$scaledPseudoD <- blocksStats$invBlockFraction * blocksStats$pseudoD
  blocksStats$scaledPseudoFd_1DD4 <- blocksStats$invBlockFraction * blocksStats$pseudoFd_1DD4
  blocksStats$scaledPseudoFd_D2D4 <- blocksStats$invBlockFraction * blocksStats$pseudoFd_D2D4
  
  # Sum up the scaled Pseudo values.
  scaledPseudo_D_Sum <- sum(blocksStats$scaledPseudoD)
  scaledPseudo_Fd_1DD4_Sum <- sum(blocksStats$scaledPseudoFd_1DD4)
  scaledPseudo_Fd_D2D4_Sum <- sum(blocksStats$scaledPseudoFd_D2D4)
  
  # Calculate the product of the global estimates and the number of blocks.
  prodGlobN_D <- nBlocks * globalEstimate_D
  prodGlobN_Fd_1DD4 <- nBlocks * globalEstimate_Fd_1DD4
  prodGlobN_Fd_D2D4 <- nBlocks * globalEstimate_Fd_D2D4
  
  fracReciprocal <- 1 / blocksStats$blockFraction - 1
  
  # Set the observed global estimate s of D, and Fd in the results object.
  fttRecord$Observed_D <- globalEstimate_D
  fttRecord$Observed_Fd_1DD4 <- globalEstimate_Fd_1DD4
  fttRecord$Observed_Fd_D2D4 <- globalEstimate_Fd_D2D4
  
  message("\t- Calculating jack-knife corrected global estimates of\nFour Taxon Test statistics.")
  # Work out the jackknife corrected estimates of D and Fds.
  fttRecord$D_jEstimate <- prodGlobN_D - scaledPseudo_D_Sum
  fttRecord$Fd_1DD4_jEstimate <- prodGlobN_Fd_1DD4 - scaledPseudo_Fd_1DD4_Sum
  fttRecord$Fd_D2D4_jEstimate <- prodGlobN_Fd_D2D4 - scaledPseudo_Fd_D2D4_Sum
  
  fttRecord$D_jVariance <- 1/nBlocks * 
    sum(1 / fracReciprocal * (((1 / blocksStats$blockFraction) * globalEstimate_D) - 
                                (fracReciprocal * blocksStats$pseudoD) - (nBlocks * globalEstimate_D) + scaledPseudo_D_Sum)^2)
  fttRecord$Fd_1DD4_jVariance <- 1/nBlocks * 
    sum(1 / fracReciprocal * (((1 / blocksStats$blockFraction) * globalEstimate_Fd_1DD4) - 
                                (fracReciprocal * blocksStats$pseudoFd_1DD4) - (nBlocks * globalEstimate_Fd_1DD4) + scaledPseudo_Fd_1DD4_Sum)^2)
  fttRecord$Fd_D2D4_jVariance <- 1/nBlocks * 
    sum(1 / fracReciprocal * (((1 / blocksStats$blockFraction) * globalEstimate_Fd_D2D4) - 
                                (fracReciprocal * blocksStats$pseudoFd_D2D4) - (nBlocks * globalEstimate_Fd_D2D4) + scaledPseudo_Fd_D2D4_Sum)^2)
  fttRecord$D_jSD <- sqrt(fttRecord$D_jVariance)
  fttRecord$Fd_1DD4_jSD <- sqrt(fttRecord$D_jVariance)
  fttRecord$Fd_D2D4_jSD <- sqrt(fttRecord$D_jVariance)
  fttRecord$D_jZ <- fttRecord$D_jEstimate / fttRecord$D_jSD
  fttRecord$Fd_1DD4_jZ <- fttRecord$Fd_1DD4_jEstimate / fttRecord$Fd_1DD4_jSD
  fttRecord$Fd_D2D4_jZ <- fttRecord$Fd_D2D4_jEstimate / fttRecord$Fd_D2D4_jSD
  fttRecord$table <- cbind(results, blocksStats)
  fttRecord$globalX2 <- -2 * sum(log(fttRecord$table$numBinomialP))
  fttRecord$X2_P <- pchisq(fttRecord$globalX2,
                           df = 2 * length(fttRecord$table$numBinomialP), 
                           lower.tail = FALSE)
  fttRecord$ABBA <- sum(blocksStats$ABBA)
  fttRecord$BABA <- sum(blocksStats$BABA)
  fttRecord$ABBAcount <- sum(blocksStats$numABBA)
  fttRecord$BABAcount <- sum(blocksStats$numBABA)
  if(fttRecord$Observed_Fd_1DD4 < 0 || is.na(fttRecord$Observed_Fd_1DD4)){fttRecord$Observed_Fd_1DD4 <- 0}
  if(fttRecord$Observed_Fd_D2D4 < 0 || is.na(fttRecord$Observed_Fd_D2D4)){fttRecord$Observed_Fd_D2D4 <- 0}
  if(fttRecord$Fd_1DD4_jEstimate < 0 || is.na(fttRecord$Fd_1DD4_jEstimate)){fttRecord$Fd_1DD4_jEstimate <- 0}
  if(fttRecord$Fd_D2D4_jEstimate < 0 || is.na(fttRecord$Fd_D2D4_jEstimate)){fttRecord$Fd_D2D4_jEstimate <- 0}
}












### BLOCK DETECTION FUNCTIONS ###

autodetect.thresholds <- function(ssdata, settings){
  ssDataTable <- ssdata$Table
  densities <- list(density(ssDataTable$AB), density(ssDataTable$AC), density(ssDataTable$BC))
  means <- list(mean(ssDataTable$AB), mean(ssDataTable$AC), mean(ssDataTable$BC))
  sds <- list(sd(ssDataTable$AB), sd(ssDataTable$AC), sd(ssDataTable$BC))
  peaks <- lapply(densities, function(x) cbind(x$x[which(diff(sign(diff(x$y))) == -2)], x$y[which(diff(sign(diff(x$y))) == -2)]))
  suitablePeaks <- lapply(1:3, function(i) which(peaks[[i]][,1] >= means[[i]] + (sds[[i]] / settings$SDstringency)))
  overallMean <- mean(c(ssDataTable$AB, ssDataTable$AC, ssDataTable$BC))
  overallSd <- sd(c(ssDataTable$AB, ssDataTable$AC, ssDataTable$BC))
  thresholds <- lapply(1:3, function(i) locate.dips(suitablePeaks[[i]], peaks[[i]], densities[[i]], settings, overallMean))
  return(thresholds)
}

# Internal function for finding interesting dips
locate.dips <- function(prospectivePeaks, peaks, density, settings, overallMean){
  if(is.integer(prospectivePeaks) && length(prospectivePeaks) > 0){
    starts <- matrix(peaks[prospectivePeaks, ], ncol=2)
    thresholds <- numeric(length=nrow(starts))
    for(i in 1:nrow(starts)){
      stopFlag <- F                                   # Set the stop flag
      last <- which(density$x == starts[i,1])         # Set the variable last as which density x is equal to the first interesting peak.
      current <- last-1                               # Set variable THIS to be one less than the last...
      while(!stopFlag){                               # As long as the stop flag is not set to true...
        if(density$y[current] < density$y[last]){     # The THIS is LESS than LAST...
          last <- current                             # Set last to this
          current <- current-1                        # Move this down one.
          if(current == 0){
            stopFlag <- T
          }
        } else {
          stopFlag <- T
        }  
      }
      thresholds[i] <- floor(density$x[last])
    }
    thresholds <- thresholds[which(thresholds > overallMean)]
    if(length(thresholds) < 1){
      if(settings$ManualFallback){
        message("\t\t\t- Falling back to manual thresholds...")
        thresholds <- settings$ManualThresholds
      } 
    }
  } else {
    if(settings$ManualFallback){
      message("\t\t\t- Falling back to manual thresholds...")
      thresholds <- settings$ManualThresholds
    } else {
      message(" - Interesting peaks was of length zero, and manual falback is off, couldn't determine thresholds")
      thresholds <- numeric(length=0)
    }
  }
  return(thresholds)
}

# Function for finding which windows meet each level of threshold.
winMeetThresh <- function(winData, thresholds){
  out <- lapply(1:length(thresholds), rep(F, times = nrow(winData)))
  if(is.numeric(thresholds)){
    out[[1]][which(winData[, 7] > thresholds[1])] <- T
    if(length(thresholds) >= 2){
      for(i in 2:length(thresholds)){
        out[[i]][which((winData[, 7] > thresholds[i]) == (winData[, 7] < thresholds[i-1]))] <- T
      }
    }
  }
  names(out) <- thresholds
  return(out)
}

runsAndLengths <- function(wmt){
  listOfruns <- lapply(wmt, function(x) rle(x))
  iter <- 1:length(listOfRuns)
  listOfBlockInd <- lapply(listOfRuns, function(x) which(x$values == T)) 
  runLengths <- lapply(iter, function(i) listOfRuns[[i]]$lengths[listOfBlockInd[[i]]])
  runLengthCumSum <- lapply(iter, function(i) cumsum(listOfRuns[[i]]$lengths)[listOfBlockInd[[i]]])
  return(list(runLens = runLengths, runCumSum = runLengthCumSum, iter = iter))
}

createBlocksDataFrames <- function(runData, scanData){
  Frames <- lapply(runData$iter, 
                     function(i){
                       data.frame(Length = runData$runLen[[i]], 
                                  Last = runData$runCumSum[[i]])
                     })
  Frames <- lapply(Frames, function(x) within(x, First <- (x$Last - (x$Length - 1))))
  # We define the start of a block as the central bp position of the first
  #  window which covers it with sufficient SS to meet the threshold.
  Frames <- lapply(Frames, function(x) within(x, FirstBP <- scanData[x$First, 4]))
  # We define the end of a block as the central bp position of the last window
  #  that covers it with a high enough SS to meet the threshold.
  Frames <- lapply(Frames, function(x) within(x, LastBP <- scanData[x$Last, 4]))
  Frames <- lapply(Frames, function(x) within(x, ApproxBpLength <- (x$LastBP - x$FirstBP) + 1))
  Frames <- lapply(Frames, function(x) x[which(x$ApproxBpLength > 1 ), ])
  Frames <- suppressWarnings(lapply(Frames, function(x) within(x, SNPs <- NA)))
  Frames <- suppressWarnings(lapply(Frames, function(x) within(x, CorrectedSNPs <- NA)))
  Frames <- suppressWarnings(lapply(Frames, function(x) within(x, P_Value <- NA)))
  Frames <- suppressWarnings(lapply(Frames, function(x) within(x, P_Threshold <- NA)))
  Frames <- suppressWarnings(lapply(Frames, function(x) within(x, fiveAge <- NA)))
  Frames <- suppressWarnings(lapply(Frames, function(x) within(x, fiftyAge <- NA)))
  Frames <- suppressWarnings(lapply(Frames, function(x) within(x, ninetyFiveAge <- NA)))
  return(Frames)
}

block.find <- function(dist, thresh){
  # Reverse the order of thresholds.
  threshRev <- rev(thresh)  
  if(is.numeric(threshRev)){
    # For each threshold, find the windows that meet it.
    windowsMeetingThresh <- winMeetThresh(dist, threshRev)
    runsData <- runsAndLengths(windowsMeetingThresh)
    blockFrames <- createBlocksDataFrames(runsData, dist)
    names(blockFrames) <- threshRev
  } else {
    Frames <- NULL
  }
  return(Frames)
}


### BLOCK DATING FUNCTIONS ###

binomcalc <- function(p, p0, N, B){pbinom(B, N, p) - p0}

#' Calculate the significance of identified recombinant blocks, and calculate their coalescence times.
#' 
#' @param blocksobj A data.frame containing all the information about the location of the identified block.
#' @param dnaobj A preprocessed DNAStringSet object.
#' @param mut A predefined mutation rate.
#' @param pthresh An alpha value threshold to use in significance testing.
#' @param bonfcorrect A logical value - whether or not to apply bonferroni correction to the pthresh.
#' @param danyway A logical value - if TRUE will esimate coalescence times of recombinant blocks even if they are not significant.
#' @param model A character value; denoting one of the mutational models allowed for in dist.dna of the APE package. 
date.blocks <- function(blocksobj, dnaobj, mut, pthresh, bonfcorrect, danyway, model){
  nBlocksToProc <- nrow(blocksobj) 
  if(!is.character(blocksobj) && nBlocksToProc > 0){
    wholeSequenceDist <- stringDist(dnaobj,
                                    method = "hamming")[1] / unique(width(dnaobj))
    if(bonfcorrect == TRUE){
      blocksobj$P_Threshold <- pthresh <- pthresh/nBlocksToProc
    } else {
      blocksobj$P_Threshold <- pthresh
    }
    for(i in 1:nBlocksToProc){
      seqBlock <- as.DNAbin(
        subseq(dnaobj, 
               start = blocksobj[i, "FirstBP"],
               end = blocksobj[i, "LastBP"])
      )
      blocksobj[i, "SNPs"] <- round(dist.dna(seqBlock, model = "N")[1])
      distanceByModel <- dist.dna(seqBlock, model = model)[1]
      blocksobj[i, "CorrectedSNPs"] <- round(distanceByModel * blocksobj[i, "ApproxBpLength"])
    }
    blocksobj$P_Value <- pbinom(blocksobj$SNPs, blocksobj$ApproxBpLength, wholeSequenceDist)
    if(!danyway){
      blocksobj <- blocksobj[which(blocksobj$P_Value < pthresh),]
      nBlocksToProc <- nrow(blocksobj)
    }
    message("\t- Checking for bad blocks. Any found will not be dated.")
    snpsBiggerThanLength <- blocksobj$ApproxBpLength <= blocksobj$CorrectedSNPs
    pValBiggerThanOne <- blocksobj$P_Value > 1
    snpsIsNA <- is.na(blocksobj$CorrectedSNPs)
    badBlocks <- which(snpsBiggerThanLength | pValBiggerThanOne | snpsIsNA)
    blocksToDate <- 1:nBlocksToProc
    if(length(badBlocks) > 0){
      blocksToDate <- blocksToDate[-badBlocks]
    }
    if(length(blocksToDate) > 0){
      for(i in blocksToDate){
        nSNPs <- blocksobj[i, "CorrectedSNPs"]
        bpLen <- blocksobj[i, "ApproxBpLength"]
        message("\t- Dating block of length: ", bpLen, " with ", nSNPs," SNPs")
        soln5 <- uniroot(binomcalc, c(0,1), p0 = 0.05,
                         B = nSNPs,
                         N = bpLen)
        soln50 <- uniroot(binomcalc, c(0,1), p0 = 0.5,
                          B = nSNPs,
                          N = bpLen)
        soln95 <- uniroot(binomcalc, c(0,1), p0 = 0.95,
                          B = nSNPs,
                          N = bpLen)
        blocksobj[i, "fiveAge"] <- round(soln5[["root"]] / (2 * mut))
        blocksobj[i, "fiftyAge"] <- round(soln50[["root"]] / (2 * mut))
        blocksobj[i, "ninetyFiveAge"] <- round(soln95[["root"]] / (2 * mut))
      }
    }
  }
  return(blocksobj)
}