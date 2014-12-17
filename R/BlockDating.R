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
                                         MutationRate <<- 10e-09
                                         PValue <<- 0.005
                                         BonfCorrection <<- TRUE
                                         DateAnyway <<- FALSE
                                         MutationCorrection <<- "JC69"
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
                                         if(!any(model == c("raw", "TS", "TV", "JC69", "K80", "F81",
                                                           "K81", "F84", "BH87", "T92", "TN93", "GG95"))){
                                           stop(paste0("Provided model must be one of the following: ", paste(c("raw", "TS", "TV", "JC69", "K80", "F81",
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
                                       },
                                     
                                     textSummary =
                                       function(){
                                         return(paste0('Settings for testing and dating recombination blocks:\n',
                                                      '-----------------------------------------------------\n',
                                                      'Assumed substitution rate for dating (MutationRate): ',
                                                      MutationRate,
                                                      '\n\nCritical alpha value for significance testing (PValue): ',
                                                      PValue,
                                                      '\n\nApply bonferroni correction to critical alpha (BonfCorrection): ',
                                                      BonfCorrection,
                                                      '\n\nKeep and date blocks that fail the alpha (DateAnyway): ', DateAnyway,
                                                      '\n\nAssumed mutation model for dating (MutationCorrection): ', MutationCorrection))
                                       },
                                     
                                     show =
                                       function(){
                                         cat(textSummary())
                                       }
                                     ))

binomcalc <- function(p, p0, N, B){pbinom(B,N,p)-p0}

date.blocks <- function(blocksobj, dnaobj, mut, pair, pthresh, bonfcorrect, danyway, model){
  if(!is.character(blocksobj) && nrow(blocksobj) > 0){ # Checking the blocksobj is not a string of characters.
    wholeSequenceDist <- stringDist(dnaobj$FullSequence[pair], method = "hamming")[1] / dnaobj$getFullLength()
    if(bonfcorrect == TRUE){
      blocksobj$P_Threshold <- pthresh <- pthresh/nrow(blocksobj)
    } else {
      blocksobj$P_Threshold <- pthresh
    }
    for(i in 1:nrow(blocksobj)){
      #Extract the two sequences required...
      Seq <- subseq(test$DNA$FullSequence[pair], start = blocksobj[i, "FirstBP"], end = blocksobj[i, "LastBP"])
      blocksobj[i, "SNPs"] <- stringDist(Seq, method = "hamming")[1]
      distanceByModel <- dist.dna(as.DNAbin(Seq), model = model)[1]
      blocksobj[i, "CorrectedSNPs"] <- round(distanceByModel * blocksobj[i, "ApproxBpLength"])
    }
    blocksobj$P_Value <- pbinom(blocksobj$SNPs, blocksobj$ApproxBpLength, wholeSequenceDist)
    if(!danyway){
      blocksobj <- blocksobj[which(blocksobj$P_Value < pthresh),]
    }
    ObservedRatio <- blocksobj$SNPs / blocksobj$ApproxBpLength
    blocksobj <- blocksobj[(blocksobj$ApproxBpLength != blocksobj$SNPs) & (blocksobj$P_Value < 1),]
    if(nrow(blocksobj) > 0){
      for(i in 1:nrow(blocksobj)){
        soln5 <- uniroot(binomcalc, c(0,1), p0 = 0.05,
                         B = blocksobj[i, "CorrectedSNPs"],
                         N = blocksobj[i, "ApproxBpLength"])
        soln50 <- uniroot(binomcalc, c(0,1), p0 = 0.5,
                          B = blocksobj[i, "CorrectedSNPs"],
                          N = blocksobj[i, "ApproxBpLength"])
        soln95 <- uniroot(binomcalc, c(0,1), p0 = 0.95,
                          B = blocksobj[i, "CorrectedSNPs"],
                          N = blocksobj[i, "ApproxBpLength"])
        blocksobj[i, "fiveAge"] <- round(soln5[["root"]] / (2 * mut))
        blocksobj[i, "fiftyAge"] <- round(soln50[["root"]]/ (2 * mut))
        blocksobj[i, "ninetyFiveAge"] <- round(soln95[["root"]] / (2 * mut))
      }
    }
  }
  return(blocksobj)
}
