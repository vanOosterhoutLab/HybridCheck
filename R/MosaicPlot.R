# Internal functions for creating the rainbow bars in ggplot2.

#' Internal function to create the vertical bars.
#' @export
vertbar_create <- function( ssobj, plotframerow, whichcomp )
{
  bool1 <- as.numeric( ssobj[,5] ) <= as.numeric( plotframerow[2] )
  bool2 <- as.numeric( plotframerow[1] ) <= as.numeric( ssobj[,6] )
  index <- which( bool1 == bool2 )
  return( mean( ssobj[index, whichcomp] ) )
}

#' Internal function to determine colours for the bars.
#' @export
col_deter <- function( invalues, reference ) 
{
  if(any( is.nan( as.numeric(invalues) ) ) ){
    warning("\nNot a numbers (NaNs)! have been detected in the plotting frame.\n
The common cause of this is a small alignment or few informative sites in the data, 
with a too high MosaicScale parameter.\nThis usually happens at either end of the 
bars and the NaNs will be dealt with my filling them in black.\n\nTo get rid of them use a lower MosaicScale parameter.")
    cols <- "#000000"
  } else {
    cols <- reference[reference[, 1] == as.numeric( invalues[1] ) & reference[,2] == as.numeric( invalues[2] ), 3]
  }
  return( cols )
}