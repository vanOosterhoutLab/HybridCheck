#' A reference class storing the settings for recombination block dating. 
#' @name BlockDatingSettings
#' @field MutationRate Numeric vector of length one. Stores the mutation rate to be used when dating blocks.
#' @field PValue Numeric vector of length one. Stores the critical alpha value for testing the signifcance of recombination regions.
#' @field BonfCorrection Logical vector of length one, stores the option of whether the critical value stored in PValue will be corrected.
#' @field DateAnyway Logical vector of length one, sotres the option of whether blocks will be dated despite failing the critical alpha.
#' @field MutationCorrection Character vector of length 1, can be any of the model supported by the ape package. Default is "HybRIDS". 
BlockDatingSettings <- setRefClass("BlockDatingSettings",
                                   fields = list(
                                     MutationRate = "numeric",
                                     PValue = "numeric",
                                     BonfCorrection = "logical",
                                     DateAnyway = "logical",
                                     MutationCorrection = "character"
                                     ),
                                   methods = list(
                                     initialize =
                                       function(){
                                         MutationRate <<- 10e-08
                                         PValue <<- 0.005
                                         BonfCorrection <<- TRUE
                                         DateAnyway <<- FALSE
                                         MutationCorrection <<- "HybRIDS"
                                       },
                                     
                                     setMutationRate =
                                       function(newRate){
                                         "Sets a new mutation rate for the settings."
                                         if(length(newRate) > 1){stop("Input must be a single value.")}
                                         MutationRate <<- newRate
                                       },
                                     
                                     setPValue =
                                       function(newValue){
                                         "Set a new critical value for significance testing of recombinant blocks."
                                         if(length(newValue) > 1){stop("Input must be a single value.")}
                                         PValue <<- newValue
                                       },
                                     
                                     setBonferonni =
                                       function(newBonf){
                                         "Set whether the critical value should be subject to bonferroni correction during block significance testing."
                                         if(length(newBonf) > 1){stop("Input must be a single value.")}
                                         BonfCorrection <<- newBonf
                                       },
                                     
                                     setDateAnyway =
                                       function(newValue){
                                         "Set whether blocks should be kept and dated even if they fail the significance test."
                                         if(length(newValue) > 1){stop("Input must be a single value.")}
                                         DateAnyway <<- newValue
                                       },
                                     
                                     setMutationCorrection =
                                       function(model){
                                         "Set the model of sequence evolution to correct the distances/number of mutations used in block dating algorithm."
                                         if(length(model) > 1){stop("Input must be a single value.")}
                                         if(!any(model == c("HybRIDS", "raw", "TS", "TV", "JC69", "K80", "F81",
                                                           "K81", "F84", "BH87", "T92", "TN93", "GG95"))){
                                           stop(paste0("Provided model must be one of the following: ", paste(c("HybRIDS", "raw", "TS", "TV", "JC69", "K80", "F81",
                                                                                                                "K81", "F84", "BH87", "T92", "TN93", "GG95."), collapse=", ")))
                                         }
                                         MutationCorrection <<- model
                                       }
                                     
                                     ))










## Internal functions for the block dating and significant values.

mergeBandD <- function(block, date) {
  output <- lapply(1:3, function(i) combineDatesWBlocks(block[[i]], date[[i]]))
  names(output) <- names(block)
  return(output)
}

combineDatesWBlocks <- function(Bs,Ds){
  outdfs <- lapply(1:length(Bs), function(i) ifelse(is.character(Bs[[i]]) || nrow(Bs[[i]]) < 1, "NO BLOCKS DETECTED OR DATED", return(comb(Bs[[i]], Ds[[i]]))))
  names(outdfs) <- names(Bs)
  return(outdfs)
}

comb <- function(B,D){
  B$fiveAge <- D[,1]
  B$fiftyAge <- D[,2]
  B$ninetyfiveAge <- D[,3]
  B$SNPnum <- D[,5]
  B$PValue <- D[,6] # Incorporate the p-value into output
  B$PThresh <- D[,7]
  B$MeanAge <- D[,"AgeMean"]
  B$CorrectedSNPs <- D[,"SNPs_Corrected"]
  B <- as.data.frame(na.omit(B))
  if(nrow(B) == 0) {
    B <- "NO BLOCKS DETECTED OR DATED"
  }
  return(B)
}

binomcalc <- function(p, p0, N, B){pbinom(B,N,p)-p0}

date.blocks <- function(blocksobj, dnaobj, mut, pair, pthresh, bonfcorrect, danyway, model){
  # Check there are blocks to date!
  if(!is.character(blocksobj)){ # Checking the blocksobj is not a string of characters.
    blocksobj <- as.matrix(blocksobj)
    blockAges <- matrix(nrow=nrow(blocksobj),ncol=9)
    wholeSequenceDist <- (dist.dna(as.DNAbin(dnaobj$FullSequence[pair,]), model="raw")[1])
    colnames(blockAges) <- c("5%","50%","95%","BlockSize","SNPs","p-value", "p-threshold", "SNPs_Corrected", "AgeMean")
    blockAges[,4] <- blocksobj[,"ApproxBpLength"]
    if( bonfcorrect == TRUE ){
      blockAges[,7] <- pthresh/nrow(blocksobj)
    } else {
      blockAges[,7] <- pthresh
    }
    for(i in 1:nrow(blocksobj)){
      #Extract the two sequences required...
      Seq <- dnaobj$FullSequence[pair,c(which(dnaobj$getFullBp() == blocksobj[i,"FirstBP"]):which(dnaobj$getFullBp() == blocksobj[i,"LastBP"]))]
      blockAges[i,5] <- dist.dna(as.DNAbin(Seq), model="N")[1]
      if(model!="HybRIDS"){
        distanceByModel <- dist.dna(as.DNAbin(Seq), model=model)[1]
        blockAges[,8] <- round(distanceByModel * blockAges[,4])
      }
    }
    blockAges[,6] <- pbinom(blockAges[,"SNPs"], blockAges[,"BlockSize"], wholeSequenceDist)
    ObservedRatio <- blockAges[,"SNPs"] / blockAges[,"BlockSize"]
    if(model == "HybRIDS"){
      ObservedRatio <- blockAges[,"SNPs"] / blockAges[,"BlockSize"]
      ActualRatio <- seq(from=0, to=2, by=0.0001)
      OutputRatio <- (0.992582633 * ActualRatio) - (0.605566567 * (ActualRatio^2) ) + (0.166571989 * (ActualRatio^3) )
      ActualRatio <- unlist(lapply(ObservedRatio, function(x) ActualRatio[which(OutputRatio >= x)][1]))
      blockAges[,8] <- round(ActualRatio * blockAges[,4])
    }
    for(i in 1:nrow(blocksobj)){
      if( blockAges[i,"p-value"] > pthresh && danyway == FALSE ) {
        next # If the block does not meet the p-value threshold, drop it and proceed to next loop iteration.
      } else {
        if(blockAges[i,"BlockSize"] == blockAges[,"SNPs"] || blockAges[i,"p-value"] == 1) next
        soln5 <- uniroot(binomcalc, c(0,1), p0=0.05, B=blockAges[i,"SNPs_Corrected"], N=blockAges[i,"BlockSize"])
        soln50 <- uniroot(binomcalc, c(0,1), p0=0.5, B=blockAges[i,"SNPs_Corrected"], N=blockAges[i,"BlockSize"])
        soln95 <- uniroot(binomcalc, c(0,1), p0=0.95, B=blockAges[i,"SNPs_Corrected"], N=blockAges[i,"BlockSize"])
        blockAges[i,3] <- round(soln5[["root"]]/(2*mut))
        blockAges[i,2] <- round(soln50[["root"]]/(2*mut))
        blockAges[i,1] <- round(soln95[["root"]]/(2*mut))
      }
    }
    blockAges[,"AgeMean"] <- ActualRatio / (2*mut) 
  } else {
    blockAges <- "NO BLOCKS TO DATE"
  }
  return(blockAges)
}
