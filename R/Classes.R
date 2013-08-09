# S3 classes for HybRIDS.
# Last Edited by B. J. Ward on 21/04/2013.

# Most of these are internal and not really for the use of the user even though some are exported... 
# Methods for them that convert dnabin classes and such from packages such as ape will be added.


# Classes defined are S3 classes.


as.HybRIDSdna <- function(x) {
  class(x) <- c("HybRIDSdna","list")
  return(x)
}



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


as.HybRIDSdatedBlocks <- function(x) {
  class(x) <- c("HybRIDSdatedBlocks","list")
  return(x)
}


as.HybRIDSdatedBlocksSET <- function(x) {
  class(x) <- c("HybRIDSdatedBlocksSET","list")
  return(x)
}

as.HybRIDSblockTable <- function(x){
  class(x) <- c("HybRIDSblockTable","data.frame")
  return(x)
}

as.HybRIDSblockTableSET <- function(x){
  class(x) <- c("HybRIDSblockTableSET","list")
  return(x)
}

as.HybRIDSdatedBlocksTable <- function(x){
  class(x) <- c("HybRIDSdatedBlocksTable","data.frame")
  return(x)
}

as.HybRIDSdatedBlocksTableSET <- function(x){
  class(x) <- c("HybRIDSdatedBlocksTableSET","list")
  return(x)
}