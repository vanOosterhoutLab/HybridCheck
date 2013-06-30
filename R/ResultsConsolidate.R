# REMOVE FROM FINAL CODE FOR DEVELOPMENT PURPOSES ONLY
load("~/Desktop/Testenv.RData")

# First we will initialize some needed variables

Consolidate.Results <- function(BlocksInput, DatesInput, intrested.in = NULL){
  unlist(lapply(BlocksInput, function(x) sum(nrow(x[[2]][[1]]), nrow(x[[2]][[1]]), nrow(x[[2]][[1]]))))
  DataFrame <- data.frame(matrix(ncol = 8))
  names(DataFrame) <- c("Pair","Beginning","End","Length","5Date","50Date","95Date","NOb")
  if ("HybRIDSblock" %in% class(BlocksInput) && "HybRIDSdates" %in% class(DatesInput)){
    cat("There isn't really any need to consolidate results for a single triplet - but let's do it anyway...")
    DataFrame <- data.frame(matrix(ncol = 8, nrow = sum(nrow(BlocksInput[[2]][[1]][[1]]), nrow(BlocksInput[[2]][[2]][[1]]), nrow(BlocksInput[[2]][[3]][[1]]))))
    names(DataFrame) <- c("Pair","Beginning","End","Length","5Date","50Date","95Date","NOb")
    DataFrame$Pair <- c(rep(names(BlocksInput[[1]][1]), nrow(BlocksInput[[2]][[1]][[1]])), rep(names(BlocksInput[[1]][2]), nrow(BlocksInput[[2]][[2]][[1]])), rep(names(BlocksInput[[1]][3]), nrow(BlocksInput[[2]][[3]][[1]])))
    DataFrame$Beginning <- c(BlocksInput[[2]][[1]][[1]]$FirstBP, BlocksInput[[2]][[2]][[1]]$FirstBP, BlocksInput[[2]][[3]][[1]]$FirstBP)
    DataFrame$End <- c(BlocksInput[[2]][[1]][[1]]$LastBP, BlocksInput[[2]][[2]][[1]]$LastBP, BlocksInput[[2]][[3]][[1]]$LastBP)
    DataFrame$Length <- c(BlocksInput[[2]][[1]][[1]]$ApproxBpLength, BlocksInput[[2]][[2]][[1]]$ApproxBpLength, BlocksInput[[2]][[3]][[1]]$ApproxBpLength)
    DataFrame[,5] <- c(unlist(lapply(DatesInput[[1]], function(x) x[,1])), unlist(lapply(DatesInput[[2]], function(x) x[,1])), unlist(lapply(DatesInput[[3]], function(x) x[,1])))
    DataFrame[,6] <- c(unlist(lapply(DatesInput[[1]], function(x) x[,2])), unlist(lapply(DatesInput[[2]], function(x) x[,2])), unlist(lapply(DatesInput[[3]], function(x) x[,2]))) 
    DataFrame[,7] <- c(unlist(lapply(DatesInput[[1]], function(x) x[,3])), unlist(lapply(DatesInput[[2]], function(x) x[,3])), unlist(lapply(DatesInput[[3]], function(x) x[,3])))  
      
    for (i in 1:3){
      
    }
    
  }
  
  
  for (i in 1:length(BlocksInput)){
    Blocks <- BlocksInput[[i]]
    Dates <- DatesInput[[i]]
    
  }
  
  
}


consolidate <- function(B, D){
  TmpDF <- data.frame(matrix(ncol = 8))
  names(TmpDF) <- c("Pair","Beginning","End","Length","5Date","50Date","95Date","NOb")
  for(i in 1:3){
    Pair <- names(B[[1]][i])
    Beginning <- B[[2]][[i]][,4]
    for (n in 1:nrow(B[[2]][i]))
    
  }
}





















