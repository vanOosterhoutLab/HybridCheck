
seq.similarity <- function(dnain, triplet, win.size, s.size, fulllength, cutbp) {
  message("Preparing input DNA sequences...")
  colnames(dnain) <- cutbp # Ensures that the correct base labels are kept.
  cutDNA <- dnain[, sequenceChecker_cpp(dnain)]
  dnaCutLength <- ncol(cutDNA)     # Measure number of bases in cut data.
  if(dnaCutLength >= 1){
    message("Checking the sliding window parameters")
    # First make sure a safe set of sliding windows is calculated for the analysis.
    # 1). Make sure the sliding window size is not some stupid value that will mess up the program loops...
    # Make sure that the window.size value entered by the user is a whole number integer.
    # If not, make it an integer.
    if(is.integer(win.size) == F){          
      win.size <- as.integer(win.size)
    }
    # Make sure the sliding window size is not bigger than the length of the dna sequence.
    # If it is, then automatically set it to 10% of the sequence length.
    # Then if setting the window size to 10% of the sequence size results in an weird window,
    # for example of size 0 or negative correct this by setting the window size to 1.
    if(win.size > dnaCutLength){
      win.size <- as.integer((dnaCutLength / 100) * 10)
      message("The set sliding window size is bigger than the length of the actual informative sites of the contig!")
      message("Default behaviour in this case is to set the sliding window to 10%
            of the sequence length... ")
      message("This is equal to ", win.size)
      triplet$SSWarning <- paste( "Sliding Window Size Warning #1: Window Size was set to", win.size, ".")
      
      if(win.size < 1L){
        win.size <- 1L
        message("The set sliding window size is bigger than the length of the actual informative sites of the contig!")
        triplet$SSWarning <- "Sliding Window Size Warning #2: Window Size was set to 1\n"
        message("Default behaviour in this case is to set the sliding window to 10%
            of the sequence length, but since this value is below 1, instead HybRIDS is
            setting the sliding window length to 1...")
      }
    }
    message("Making all the window frames...")
    if(win.size >= 1L) {
      allstepsfrom <- 1 + as.integer(win.size / 2)
      allstepsto <- (ncol(cutDNA) - as.integer(win.size / 2)) + 1
      allsteps <- seq(from = allstepsfrom, to = allstepsto, by = s.size)
      windowp1 <- allsteps - as.integer(win.size / 2) # All the window start points.
      windowp2 <- allsteps + as.integer(win.size / 2) # All the window end points.
      removals <- which(windowp2 > dnaCutLength) # Remove the last window and any accidentally beyond the sequence end point.
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
      Distances[, 4] <- as.numeric(unlist(lapply(1:length(allsteps), function(i) colnames(cutDNA)[allsteps[i]]))) # ActualBP Center
      Distances[, 5] <- as.numeric(colnames(cutDNA)[windowp1]) # Actual BP Start
      Distances[, 6] <- as.numeric(colnames(cutDNA)[windowp2]) # Actual BP End
      rm(windowp1, windowp2, allsteps, allstepsto, allstepsfrom)
      colnames(Distances) <- c("WindowCenter", "WindowStart", "WindowEnd", "ActualCenter", "ActualStart", "ActualEnd", unlist(lapply(pairs, function(x) paste(LETTERS[x], collapse=""))))
      # Set up the loop for calculation.
      message("Analysing Now!")
      #n1 <- 0
      #Do the loop - Calculates all the hamming distances for all contig pairs, in all window frames. 
      for(i in seq(nrow(Distances))) {
        dnaStretch <- cutDNA[, Distances[i, 2] : Distances[i, 3]]
        Distances[i, 7] <- sum(dnaStretch[1, ] != dnaStretch[2, ])
        Distances[i, 8] <- sum(dnaStretch[1, ] != dnaStretch[3, ])
        Distances[i, 9] <- sum(dnaStretch[2, ] != dnaStretch[3, ])
      }      
      Distances[ , c(7, 8, 9)] <- 100 - round((as.numeric(Distances[ , c(7, 8, 9)]) / (win.size + 1)) * 100)
      triplet$SSTable <- as.data.frame(Distances)
      # Ok Triplet SS table made successfully so let's get rid of any error messages to do with it, since everything was fine.
      triplet$SSError <- character()
      } else {
        triplet$SSError <- c(triplet$SSError,"The sliding window size is less than 1, this is not supposed to happen")
      }
      triplet$WindowSizeUsed <- win.size
      triplet$StepSizeUsed <- s.size
      triplet$InformativeDNALength <- dnaCutLength
    } else {
      triplet$SSError <- c(triplet$SSError,"There are no informative sites to work on - skipping analysis of this triplet...")
    }
}