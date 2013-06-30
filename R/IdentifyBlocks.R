# Functions for Identifying potential blocks from sequence similarity data, and dating them.
# Created by Ben J. Ward 30/04/2013.
# Last Edited by Ben J. Ward on 30/04/2013.

Identify.Blocks <- function(similarity.data, manual.thresholds=90, autodetect.thresholds=T, use.manual.fallback=T, standard.deviation.stringency=2) {
  stopifnot(is.numeric(manual.thresholds),is.logical(autodetect.thresholds),is.logical(use.manual.fallback), is.numeric(standard.deviation.stringency))
  if("HybRIDSseqsim" %in% class(similarity.data)){
    Blocks <- identify.blocks(similarity.data, manual.thresholds, autodetect.thresholds, use.manual.fallback, standard.deviation.stringency)
    return(Blocks)
  } else {
    if("HybRIDSseqsimSET" %in% class(similarity.data)){
      Blocks <- lapply(similarity.data, function(x) identify.blocks(x, manual.thresholds, autodetect.thresholds, use.manual.fallback, standard.deviation.stringency))
      return(as.HybRIDSblockSET(Blocks))
    } else {
      stop("You did not enter a valid datatype as argument similarity.data")
    }
  }
}






autodetect.thresholds <- function(data, SDdiv, manual, manoveride) {
  Densities <- lapply(7:9, function(i) density(data[,i])) # Generate the Densities for each of the three triplet comparrisons.
  Means <- lapply(7:9, function(i) mean(data[,i]))        # Generate the Densities for each of the three triplet comparrisons.
  Sds <- lapply(7:9, function(i) sd(data[,i]))            # Generate the Standard Deviations for each of the three triplet comparrisons.
  Highinds <- lapply(Densities, function(x) which(diff(sign(diff(x$y))) == -2))  # Find the Peaks of the Densities.
  Peaks <- lapply(Densities, function(x) cbind(x$x[which(diff(sign(diff(x$y))) == -2)],x$y[which(diff(sign(diff(x$y))) == -2)]))
  SuitablePeaks <- lapply(1:3, function(i) which(Peaks[[i]][,1] >= Means[[i]]+(Sds[[i]]/SDdiv)))
  Thresholds <- lapply(1:3, function(i) locate.dips(SuitablePeaks[[i]],Peaks[[i]],Densities[[i]],manual,manoveride)) # Use the SuitablePeaks identified to find the low points preceeding them to get the thresholds. Uses interesting.dips function.
  Thresholds <- lapply(Thresholds, function(x) x[which(x > mean(c(data[,7],data[,8],data[,9])))])
  Thresholds <- lapply(Thresholds, function(x) ifelse(length(x)>0,x,"NO SUITABLE THRESHOLDS RESULTED FROM POTENTIAL PEAKS"))
  return(Thresholds)
}



locate.dips <- function(interesting.peaks, peaks, density, man, manfall) {
  if(is.integer(interesting.peaks) && length(interesting.peaks) > 0){
    Starts <- matrix(peaks[interesting.peaks,],ncol=2)
    thresholds <- numeric()
    for(i in 1:nrow(Starts)){
      Stop <- F                                  # Set the stop flag
      Last <- which(density$x == Starts[i,1])    # Set the variable last as which density x is equal to the first interesting peak.
      This <- Last-1                             # Set variable THIS to be one less than the last...
      while(Stop != T){                          # As long as the stop flag is not set to true...
        if(density$y[This] < density$y[Last]){   # The THIS is LESS than LAST...
          Last <- This                           # Set last to this
          This <- This-1                         # Move this down one.
        } else {
          Stop <- T
        }  
      }
      thresholds[i] <- floor(density$x[Last])
    }
  } else {
    if(manfall==T){
      cat("Falling back to manual thresholds...\n")
      thresholds <- floor(man)
    } else {
      thresholds <- "INTERESTING PEAKS WAS OF LENGTH ZERO, COULDN'T DETERMINE THRESHOLDS"
    }
  }
  return(thresholds)
}



identify.blocks <- function(x, manual.thresholds=90, autodetect=T, manualfallback=T, SDstringency=2) {
  stopifnot(is.logical(autodetect),is.logical(manualfallback),is.numeric(SDstringency),"HybRIDSseqsim" %in% class(x), is.numeric(manual.thresholds))
  distances <- x$Distances
  if(autodetect==T) {
    cat("Using the autodetect thresholds method...\n")
    cat("Deciding on suitable thresholds...\n")
    Thresholds <- autodetect.thresholds(distances,SDstringency,manual.thresholds,manualfallback)     # Autodetect the thresholds for block identification.. uses internal function.
    names(Thresholds) <- c(paste(x$ContigNames[1],x$ContigNames[2],sep=":"),paste(x$ContigNames[1],x$ContigNames[3],sep=":"),paste(x$ContigNames[2],x$ContigNames[3],sep=":"))
  } else {
    Thresholds <- list(manual.thresholds,manual.thresholds,manual.thresholds)
    names(Thresholds) <- names(Thresholds) <- c(paste(x$ContigNames[1],x$ContigNames[2],sep=":"),paste(x$ContigNames[1],x$ContigNames[3],sep=":"),paste(x$ContigNames[2],x$ContigNames[3],sep=":"))
  }
  cat("Now beginning Block Search...\n\n")
  Blocks <- lapply(1:3, function(i) block.find2(distances[,c(1:6,6+i)], Thresholds[[i]]))
  names(Blocks) <- names(Thresholds) <- c(paste(x$ContigNames[1],x$ContigNames[2],sep=":"),paste(x$ContigNames[1],x$ContigNames[3],sep=":"),paste(x$ContigNames[2],x$ContigNames[3],sep=":"))
  
  return(as.HybRIDSblock(list(Thresholds,Blocks)))
}



block.find2 <- function(dist,thresh) {
  thresh2 <- rev(thresh)                                                           # Reverse the order of thresholds.
  b1 <- list()                                                                     # Create and empty list called b1.
  length(b1) <- length(thresh2)
  
  if(is.numeric(thresh2)){
    
    for(i in 1:length(thresh2)){                                                     #Find which sliding windows meet the threshold with for loop.
      b1[[i]] <- rep(F, times=nrow(dist))
      if(i == 1){
        b1[[i]][which(dist[,7] > thresh2[i])] <- T
      } else {
        b1[[i]][which((dist[,7] > thresh2[i]) == (dist[,7] < thresh2[i-1]))] <- T
      }
    }
    names(b1) <- thresh2
    runs <- lapply(b1, function(x) rle(x))                                           # Generate a list of runs data...
    ind <- lapply(runs, function(x) which(x$values==T))                              # Identify which runs are blocks and not runs of non-blocks.
    runs2 <- lapply(1:length(runs), function(i) runs[[i]]$lengths[ind[[i]]])         # Get the lengths of all the runs for T, meaning blocks.
    sums <- lapply(1:length(runs), function(i) cumsum(runs[[i]]$lengths)[ind[[i]]])  # Generate a cumulative sum of the run lengths and pick out the ones for the actual blocks.
    BlockPos2 <- lapply(1:length(runs), function(i) data.frame(Length=runs2[[i]], Last=sums[[i]]))
    BlockPos2 <- lapply(BlockPos2, function(x) within(x, First <- (x$Last - (x$Length-1))))
    BlockPos2 <- lapply(BlockPos2, function(x) within(x, FirstBP <- dist[x$First,5]))
    BlockPos2 <- lapply(BlockPos2, function(x) within(x, LastBP <- dist[x$Last,6]))
    BlockPos2 <- lapply(BlockPos2, function(x) within(x, ApproxBpLength <- (x$LastBP-x$FirstBP)+1))
    names(BlockPos2) <- thresh2
    
  } else {
    BlockPos2 <- "NO SUITABLE THRESHOLDS WITH WHICH TO FIND BLOCKS"
  }
  return(BlockPos2)
}