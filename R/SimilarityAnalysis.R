scan.similarity <- function(dna, triplet, settings){
  message(paste0("Scanning sequence similarity for triplet ", paste0(triplet$ContigNames, collapse=", ")))
  dnain <- dna$pullTriplet(triplet$ContigNames)
  precons <- which(colSums(consensusMatrix(dnain) != 0) > 1)
  cutDNA <- DNAStringSet(character(length = 3))
  cutDNA[[1]] <- dnain[[1]][precons]
  cutDNA[[2]] <- dnain[[2]][precons]
  cutDNA[[3]] <- dnain[[3]][precons]
  names(cutDNA) <- names(dnain)
  applicableBP <- dna$InformativeBp[precons]
  triplet$readSettings(cutDNA, settings)
  if(triplet$InformativeDNALength >= 1){
    message("Checking the sliding window parameters...")
    if(triplet$ScanData$WindowSizeUsed > triplet$InformativeDNALength){
      triplet$ScanData$WindowSizeUsed <- as.integer((triplet$InformativeDNALength / 100) * 10)
      message("The set sliding window size is bigger than the length of the actual informative sites of the contig!")
      message("Continuing with analysis but set the sliding window to 10%
            of the sequence length... ")
      message("This is equal to ", triplet$ScanData$WindowSizeUsed)
      if(triplet$ScanData$WindowSizeUsed < 1L){
        triplet$ScanData$WindowSizeUsed <- 1L
        message("Default behaviour in this case is to set the sliding window to 10%
            of the sequence length, but since this value is below 1, instead setting
                the sliding window length to 1...")
      }
    }
    message("Making all the window frames...")
    if(triplet$ScanData$WindowSizeUsed >= 1L) {
      halfWindow <- as.integer(triplet$ScanData$WindowSizeUsed / 2)
      allstepsfrom <- 1 + halfWindow
      allstepsto <- (triplet$InformativeDNALength - halfWindow) + 1
      #allsteps <- seq(from = allstepsfrom, to = allstepsto, by = triplet$ScanData$StepSizeUsed)
      #windowp1 <- allsteps - halfWindow # All the window start points.
      #windowp2 <- allsteps + halfWindow # All the window end points.
      windowp2 <- seq(from = triplet$ScanData$WindowSizeUsed, to = triplet$InformativeDNALength, by = triplet$ScanData$StepSizeUsed)
      windowp1 <- windowp2 - (triplet$ScanData$WindowSizeUsed - 1)
      allsteps <- floor((windowp1 + windowp2) / 2)
      
      removals <- which(windowp2 > triplet$InformativeDNALength) # Remove the last window and any accidentally beyond the sequence end point.
      if(length(removals) > 0) {
        allsteps <- allsteps[-removals]
        windowp1 <- windowp1[-removals]
        windowp2 <- windowp2[-removals]
      }
      pairs <- combn(1:3 , 2, simplify = F) # Generate all triplets pairs.
      Distances <- matrix(ncol = 9, nrow = length(windowp1))
      Distances[, 1] <- allsteps
      Distances[, 2] <- windowp1
      Distances[, 3] <- windowp2
      Distances[, 4] <- as.numeric(unlist(lapply(1:length(allsteps), function(i) applicableBP[allsteps[i]]))) # ActualBP Center
      Distances[, 5] <- as.numeric(applicableBP[windowp1]) # Actual BP Start
      Distances[, 6] <- as.numeric(applicableBP[windowp2]) # Actual BP End
      rm(windowp1, windowp2, allsteps, allstepsto, allstepsfrom)
      colnames(Distances) <- c("WindowCenter", "WindowStart", "WindowEnd", "ActualCenter", "ActualStart", "ActualEnd", unlist(lapply(pairs, function(x) paste(LETTERS[x], collapse=""))))
      # Set up the loop for calculation.
      message("Scanning Now!")
      conMatAB <- colSums(consensusMatrix(cutDNA[c(1, 2)]) != 0) > 1
      conMatAC <- colSums(consensusMatrix(cutDNA[c(1, 3)]) != 0) > 1
      conMatBC <- colSums(consensusMatrix(cutDNA[c(2, 3)]) != 0) > 1
      # Do the loop - Calculates all the hamming distances for all contig pairs, in all window frames.
      for(i in seq(nrow(Distances))){
        stretch <- Distances[i, 2] : Distances[i, 3]
        Distances[i, 7] <- sum(conMatAB[stretch])
        Distances[i, 8] <- sum(conMatAC[stretch])
        Distances[i, 9] <- sum(conMatBC[stretch])
      }      
      Distances[ , c(7, 8, 9)] <- 100 - round((as.numeric(Distances[ , c(7, 8, 9)]) / (triplet$ScanData$WindowSizeUsed + 1)) * 100)
      triplet$ScanData$Table <- as.data.frame(Distances)
      } else {
        stop("The sliding window size is less than 1, this is not supposed to be possible.")
      }
    } else {
      warning(paste0("There are no informative sites to work on - skipping analysis of this triplet: ", triplet$ContigNames, collapse=", "))
    }
}