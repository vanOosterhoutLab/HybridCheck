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
    cols <- "#000000"
  } else {
    cols <- reference[reference[, 1] == as.numeric( invalues[1] ) & reference[,2] == as.numeric( invalues[2] ), 3]
  }
  return( cols )
}