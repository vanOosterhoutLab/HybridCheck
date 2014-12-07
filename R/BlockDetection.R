# Functions and classes for detecting blocks.
#' Reference type for storing the settings for block detection.
#' @name BlockDetectionSettings
#' @field ManualThresholds A numeric value between 1 and 100, the percentage raw sequence similarity a region has to reach, before it is identified as a block.
BlockDetectionSettings <- setRefClass("BlockDetectionSettings",
                                      
                                      fields = list(
                                        ManualThresholds = "numeric",
                                        AutoThresholds = "logical",
                                        ManualFallback = "logical",
                                        SDstringency = "numeric"
                                        ),
                                      
                                      methods = list(
                                        initialize = function(){
                                          "Initializes the settingsobject."
                                          ManualThresholds <<- 90
                                          AutoThresholds <<- TRUE
                                          ManualFallback <<- TRUE
                                          SDstringency <<- 2
                                        },
                                        
                                        setManualThresholds =
                                          function(values){
                                            "Checks input values for changing the manual sequence similarity thresholds for block detection and sets the parmeter."
                                            if(any(values > 100) || any(values < 0)){stop("Enter a numeric value between 1 and 100.")}
                                            ManualThresholds <<- values
                                          },
                                        
                                        setSDstringency =
                                          function(value){
                                            "Checks input values for changing the SD stringency parameter for block detection and sets the parameter."
                                            if(value == 0){stop("You can't enter a zero value.")}
                                            SDstringency <<- value
                                          },
                                        
                                        setSettings =
                                          function(...){
                                            settings <- list(...)
                                            parameters <- names(settings)
                                            for(i in 1:length(settings)){
                                              if(parameters[i] == "ManualThresholds"){
                                                setManualThresholds(settings[[i]])
                                              }
                                              if(parameters[i] == "AutoThresholds"){
                                                AutoThresholds <<- settings[[i]]
                                              }
                                              if(parameters[i] == "ManualFallback"){
                                                ManualFallback <<- settings[[i]]
                                              }
                                              if(parameters[i] == "SDstringency"){
                                                setSDstringency(settings[[i]])
                                              }
                                            }
                                          },
                                        
                                        textSummary =
                                          function(){
                                            return(paste0('Settings for detecting blocks from recombination signal:\n',
                                                         '--------------------------------------------------------\n',
                                                         'Manual sequence similarity thresholds (ManualThresholds): ',
                                                         paste(ManualThresholds, collapse=", "),
                                                         '\n\nAutomatically decide on thresholds (AutoThresholds): ',
                                                         AutoThresholds,
                                                         '\n\nFall back to manual thresholds (ManualFallback): ',
                                                         ManualFallback,
                                                         '\n\nStandard deviation divisor (SDStringency): ', SDstringency
                                              ))
                                          },
                                        
                                        show =
                                          function(){
                                            cat(textSummary())
                                          }
                                        )
                                      )



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
        message("Falling back to manual thresholds...")
        thresholds <- settings$ManualThresholds
      } 
    }
  } else {
    if(settings$ManualFallback){
      message("Falling back to manual thresholds...")
      thresholds <- settings$ManualThresholds
    } else {
      message("Interesting peaks was of length zero, and manual falback is off, couldn't determine thresholds")
      thresholds <- numeric(length=0)
    }
  }
  return(thresholds)
}

block.find <- function(dist,thresh){
  thresh2 <- rev(thresh)                                                           # Reverse the order of thresholds.
  b1 <- list()                                                                     # Create and empty list called b1.
  length(b1) <- length(thresh2)
  if(is.numeric(thresh2)){
    for(i in 1:length(thresh2)){                                                   #Find which sliding windows meet the threshold with for loop.
      b1[[i]] <- rep(F, times = nrow(dist))
      if(i == 1){
        b1[[i]][which(dist[,7] > thresh2[i])] <- T
      } else {
        b1[[i]][which((dist[,7] > thresh2[i]) == (dist[,7] < thresh2[i-1]))] <- T
      }
    }
    names(b1) <- thresh2
    runs <- lapply(b1, function(x) rle(x))                                           # Generate a list of runs data...
    ind <- lapply(runs, function(x) which(x$values == T))                              # Identify which runs are blocks and not runs of non-blocks.
    runs2 <- lapply(1:length(runs), function(i) runs[[i]]$lengths[ind[[i]]])         # Get the lengths of all the runs for T, meaning blocks.
    sums <- lapply(1:length(runs), function(i) cumsum(runs[[i]]$lengths)[ind[[i]]])  # Generate a cumulative sum of the run lengths and pick out the ones for the actual blocks.
    BlockPos2 <- lapply(1:length(runs), function(i) data.frame(Length = runs2[[i]], Last = sums[[i]]))
    BlockPos2 <- lapply(BlockPos2, function(x) within(x, First <- (x$Last - (x$Length - 1))))
    BlockPos2 <- lapply(BlockPos2, function(x) within(x, FirstBP <- dist[x$First,4])) # We define the start of a block as the central bp position of the first window which covers it with sufficient SS to meet the threshold.
    BlockPos2 <- lapply(BlockPos2, function(x) within(x, LastBP <- dist[x$Last,4])) # We define the end of a block as the central bp position of the last window that covers it with a high enough SS to meet the threshold.
    BlockPos2 <- lapply(BlockPos2, function(x) within(x, ApproxBpLength <- (x$LastBP - x$FirstBP) + 1))
    BlockPos2 <- lapply(BlockPos2, function(x) x[which(x$ApproxBpLength > 1 ),]) # Remove any blocks with a bpsize of 1 straight away.
    BlockPos2 <- suppressWarnings(lapply(BlockPos2, function(x) within(x, SNPs <- NA)))
    BlockPos2 <- suppressWarnings(lapply(BlockPos2, function(x) within(x, CorrectedSNPs <- NA)))
    BlockPos2 <- suppressWarnings(lapply(BlockPos2, function(x) within(x, P_Value <- NA)))
    BlockPos2 <- suppressWarnings(lapply(BlockPos2, function(x) within(x, P_Threshold <- NA)))
    BlockPos2 <- suppressWarnings(lapply(BlockPos2, function(x) within(x, fiveAge <- NA)))
    BlockPos2 <- suppressWarnings(lapply(BlockPos2, function(x) within(x, fiftyAge <- NA)))
    BlockPos2 <- suppressWarnings(lapply(BlockPos2, function(x) within(x, ninetyFiveAge <- NA)))
    names(BlockPos2) <- thresh2
  } else {
    BlockPos2 <- NULL
  }
  return(BlockPos2)
}