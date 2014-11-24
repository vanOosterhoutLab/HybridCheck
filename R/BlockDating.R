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
                                       },
                                     
                                     setSettings =
                                       function(...){
                                         settings <- list(...)
                                         parameters <- names(settings)
                                         for(i in 1:length(settings)){
                                           if(parameters[i] == "MutationCorrection"){
                                             setMutationCorrection(settings[[i]])
                                           }
                                           if(parameters[i] == "DateAnyway"){
                                             setDateAnyway(settings[[i]])
                                           }
                                           if(parameters[i] == "BonfCorrection"){
                                             setBonferonni(settings[[i]])
                                           }
                                           if(parameters[i] == "PValue"){
                                             setPValue(settings[[i]])
                                           }
                                           if(parameters[i] == "MutationRate"){
                                             setMutationRate(settings[[i]])
                                           }
                                         }
                                       }
                                     
                                     ))

binomcalc <- function(p, p0, N, B){pbinom(B,N,p)-p0}

date.blocks <- function(blocksobj, dnaobj, mut, pair, pthresh, bonfcorrect, danyway, model){
  if(!is.character(blocksobj)){ # Checking the blocksobj is not a string of characters.
    wholeSequenceDist <- (dist.dna(as.DNAbin(dnaobj$FullSequence[pair,]), model="raw")[1])
    if(bonfcorrect == TRUE){
      blocksobj$P_Threshold <- pthresh <- pthresh/nrow(blocksobj)
    } else {
      blocksobj$P_Threshold <- pthresh
    }
    for(i in 1:nrow(blocksobj)){
      #Extract the two sequences required...
      Seq <- dnaobj$FullSequence[pair,c(which(dnaobj$getFullBp() == blocksobj[i,"FirstBP"]):which(dnaobj$getFullBp() == blocksobj[i,"LastBP"]))]
      blocksobj[i,"SNPs"] <- dist.dna(as.DNAbin(Seq), model="N")[1]
      if(model!="HybRIDS"){
        distanceByModel <- dist.dna(as.DNAbin(Seq), model=model)[1]
        blocksobj[i,"CorrectedSNPs"] <- round(distanceByModel * blocksobj[i,"ApproxBpLength"])
      }
    }
    blocksobj$P_Value <- pbinom(blocksobj$SNPs, blocksobj$ApproxBpLength, wholeSequenceDist)
    if(!danyway){
      blocksobj <- blocksobj[which(blocksobj$P_Value < pthresh),]
    }
    ObservedRatio <- blocksobj$SNPs / blocksobj$ApproxBpLength
    if(model == "HybRIDS"){
      ActualRatio <- seq(from=0, to=2, by=0.0001)
      OutputRatio <- (0.992582633 * ActualRatio) - (0.605566567 * (ActualRatio^2) ) + (0.166571989 * (ActualRatio^3))
      ActualRatio <- unlist(lapply(ObservedRatio, function(x) ActualRatio[which(OutputRatio >= x)][1]))
      blocksobj$CorrectedSNPs <- round(ActualRatio * blocksobj$ApproxBpLength)
    }
    blocksobj <- blocksobj[(blocksobj$ApproxBpLength != blocksobj$SNPs) & (blocksobj$P_Value < 1),]
    for(i in 1:nrow(blocksobj)){
        soln5 <- uniroot(binomcalc, c(0,1), p0=0.05, B=blocksobj[i,"CorrectedSNPs"], N=blocksobj[i,"ApproxBpLength"])
        soln50 <- uniroot(binomcalc, c(0,1), p0=0.5, B=blocksobj[i,"CorrectedSNPs"], N=blocksobj[i,"ApproxBpLength"])
        soln95 <- uniroot(binomcalc, c(0,1), p0=0.95, B=blocksobj[i,"CorrectedSNPs"], N=blocksobj[i,"ApproxBpLength"])
        blocksobj[i,"fiveAge"] <- round(soln5[["root"]]/(2*mut))
        blocksobj[i,"fiftyAge"] <- round(soln50[["root"]]/(2*mut))
        blocksobj[i,"ninetyFiveAge"] <- round(soln95[["root"]]/(2*mut))
    }
  }
  return(blocksobj)
}
