# Internal functions for putative block detection.

autodetect.thresholds <- function( ssdatatable, SDdiv, manual, manoveride ) {
  
  # Internal function for finding interesting dips
  locate.dips <- function( interesting.peaks, peaks, density, man, manfall, noisemean ) {
    if( is.integer( interesting.peaks ) && length( interesting.peaks ) > 0 ) {
      Starts <- matrix( peaks[interesting.peaks, ], ncol=2 )
      thresholds <- numeric( length=nrow( Starts ) )
      for( i in 1:nrow( Starts ) ) {
        Stop <- F                                     # Set the stop flag
        Last <- which( density$x == Starts[i,1] )     # Set the variable last as which density x is equal to the first interesting peak.
        This <- Last-1                                # Set variable THIS to be one less than the last...
        while( Stop != T ) {                          # As long as the stop flag is not set to true...
          if( density$y[This] < density$y[Last] ) {   # The THIS is LESS than LAST...
            Last <- This                              # Set last to this
            This <- This-1                            # Move this down one.
            if( This == 0 ) {
              Stop <- T
            }
          } else {
            Stop <- T
          }  
        }
        thresholds[i] <- floor( density$x[Last] )
      }
      thresholds <- thresholds[which( thresholds > noisemean )]
      if( length( thresholds ) < 1 ) {
        if( manfall==T ) {
          cat( "Falling back to manual thresholds...\n" )
          thresholds <- floor( man )
        } else {
          cat( "Valid thresholds could not be auto determined from suitable peaks\nfallback to manual thresholds is off" )
          thresholds <- "VALID THRESHOLDS COULD NOT BE AUTO DETERMINED FROM SUITABLE PEAKS\nFALLBACK TO MANUAL THRESHOLDS IS OFF"
        }
      }
    } else {
      if( manfall == T ) {
        cat( "Falling back to manual thresholds...\n" )
        thresholds <- floor( man )
      } else {
        cat( "Interesting peaks was of length zero, and manual falback is off, couldn't determine thresholds\n" )
        thresholds <- "INTERESTING PEAKS WAS OF LENGTH ZERO, AND MANUAL FALLBACK IS OFF COULDN'T DETERMINE THRESHOLDS"
      }
    }
    return( thresholds )
  } # END OF INTERNAL FUNCTION.
  
  Densities <- lapply( 7:9, function(i) density( ssdatatable[,i] ) ) # Generate the Densities for each of the three triplet comparrisons.
  Means <- lapply( 7:9, function(i) mean( ssdatatable[,i] ) )        # Generate the Means for each of the three triplet comparrisons.
  Sds <- lapply( 7:9, function(i) sd( ssdatatable[,i] ) )            # Generate the Standard Deviations for each of the three triplet comparrisons.
  # Find the Peaks of the Densities.
  Peaks <- lapply( Densities, function(x) cbind( x$x[ which( diff( sign( diff( x$y ) ) ) == -2 ) ], x$y[ which( diff( sign( diff( x$y ) ) ) == -2 ) ] ) )
  SuitablePeaks <- lapply( 1:3, function(i) which( Peaks[[i]][,1] >= Means[[i]] + ( Sds[[i]] / SDdiv ) ) )
  # Use of the dip location internal function here...
  noisymean <- mean( c( ssdatatable[,7], ssdatatable[,8], ssdatatable[,9] ) )
  noisysd <- sd( c( ssdatatable[,7], ssdatatable[,8], ssdatatable[,9] ) )
  Thresholds <- lapply( 1:3, function(i) locate.dips( SuitablePeaks[[i]], Peaks[[i]], Densities[[i]], manual, manoveride, noisymean) ) # Use the SuitablePeaks identified to find the low points preceeding them to get the thresholds. Uses interesting.dips function.
  return( Thresholds )
}

# Internal block detection function.
block.find <- function( dist,thresh ) {
  thresh2 <- rev( thresh )                                                           # Reverse the order of thresholds.
  b1 <- list()                                                                     # Create and empty list called b1.
  length( b1 ) <- length( thresh2 )
  if( is.numeric( thresh2 ) ) {
    for(i in 1:length( thresh2 ) ) {                                                   #Find which sliding windows meet the threshold with for loop.
      b1[[i]] <- rep( F, times = nrow( dist ) )
      if( i == 1 ) {
        b1[[i]][which( dist[,7] > thresh2[i] )] <- T
      } else {
        b1[[i]][which( ( dist[,7] > thresh2[i] ) == ( dist[,7] < thresh2[i-1] ) )] <- T
      }
    }
    names( b1 ) <- thresh2
    runs <- lapply( b1, function(x) rle( x ) )                                           # Generate a list of runs data...
    ind <- lapply( runs, function(x) which( x$values == T ) )                              # Identify which runs are blocks and not runs of non-blocks.
    runs2 <- lapply( 1:length( runs ), function(i) runs[[i]]$lengths[ind[[i]]] )         # Get the lengths of all the runs for T, meaning blocks.
    sums <- lapply( 1:length( runs ), function(i) cumsum( runs[[i]]$lengths )[ind[[i]]] )  # Generate a cumulative sum of the run lengths and pick out the ones for the actual blocks.
    BlockPos2 <- lapply( 1:length( runs ), function(i) data.frame( Length = runs2[[i]], Last = sums[[i]] ) )
    BlockPos2 <- lapply( BlockPos2, function(x) within( x, First <- ( x$Last - ( x$Length - 1 ) ) ) )
    BlockPos2 <- lapply( BlockPos2, function(x) within( x, FirstBP <- dist[x$First,4] ) ) # We define the start of a block as the central bp position of the first window which covers it with sufficient SS to meet the threshold.
    BlockPos2 <- lapply( BlockPos2, function(x) within( x, LastBP <- dist[x$Last,4] ) ) # We define the end of a block as the central bp position of the last window that covers it with a high enough SS to meet the threshold.
    BlockPos2 <- lapply( BlockPos2, function(x) within( x, ApproxBpLength <- ( x$LastBP - x$FirstBP ) + 1 ) )
    BlockPos2 <- lapply( BlockPos2, function(x) x[which( x$ApproxBpLength > 1 ),] ) # Remove any blocks with a bpsize of 1 straight away.
    names( BlockPos2 ) <- thresh2
    for( i in 1:length( BlockPos2 ) ) {
      if( nrow( BlockPos2[[i]] ) < 1) {
        BlockPos2[[i]] <- "NO BLOCKS DETECTED UNDER THIS THRESHOLD"
      }
    }
  } else {
    BlockPos2 <- "NO SUITABLE THRESHOLD TO ID BLOCKS WITH"
  }
  return( BlockPos2 )
}