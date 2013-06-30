# S3 classes for HybRIDS.
# Last Edited by B. J. Ward on 21/04/2013.




# Classes defined are S3 classes.

#' DNA sequence data class for HybRIDS.
#' 
as.HybRIDSdna <- function(x) {
  class(x) <- c("HybRIDSdna","list")
  return(x)
}


#' A results class for HybRIDS output.
#' 
as.HybRIDSseqsim <- function(x) {
  class(x) <- c("HybRIDSseqsim","list")
  return(x)
}


as.HybRIDSseqsimSET <- function(x) {
  class(x) <- c("HybRIDSseqsimSET","list")
  return(x)
}


as.HybRIDSblock <- function(x) {
  class(x) <- c("HybRIDSblock","list")
  return(x)
}


as.HybRIDSblockSET <- function(x) {
  class(x) <- c("HybRIDSblockSET","list")
  return(x)
}


as.HybRIDSlinesplot <- function(x) {
  class(x) <- c("HybRIDSlinesplot","gg","ggplot")
  return(x)
}


as.HybRIDSbars <- function(x) {
  class(x) <- c("HybRIDSbars","list")
  return(x)
}


as.HybRIDSdates <- function(x) {
  class(x) <- c("HybRIDSdates","list")
  return(x)
}


as.HybRIDSdatesSET <- function(x) {
  class(x) <- c("HybRIDSdatesSET","list")
  return(x)
}