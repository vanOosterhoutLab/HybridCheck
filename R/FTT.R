#' A Reference Class for performing and storing results of four taxon tests for introgression.
#' @name FTTester
FTTester <- setRefClass("FTTester",
                         
                         fields = list(global = "logical",
                                       jackKnife = "logical",
                                       jackKnifeBlockSize = "integer",
                                       
                                       
                                       local = "logical",
                                       windowSize = "integer",
                                       stepSize = "integer",
                                       
                                       
                                       results = "list"),
                         
                         methods = list(
                           
                           runTest =
                             function(dna, pops){
                               
                               # Calculate the number of alleles at each site in the alignment and
                               # find which ones are bi-allelic, and which ones are variable.
                               numberOfAlleles.all <- getNumberOfAlleles(dna)
                               biAllelicSites.all <- getBiAllelicSites(numberOfAlleles.all)
                               variableSites <- getVariableSites(numberOfAlleles.all)
                               # Make a version of the sequence alignment which only includes variable sites,
                               # and get the number of alleles at those variable sites,
                               # find which of those variable sites are bi-allelic.
                               variableAlignment <- subsetSequence(dna, variableSites)
                               numberOfAlleles.var <- numberOfAlleles.all[variableSites]
                               S.all <- length(numberOfAlleles.var)
                               biAllelicSites.var <- which(numberOfAlleles.var == 2)
                               
                               # Population slice calculates for a subset of the alignment of variable sites:
                               # The counts, the number of alleles, and the counts of the (global)bi-allelic sites. 
                               p1Slice <- populationSlice(variableAlignment[pops[[1]]], biAllelicSites.var)
                               p2Slice <- populationSlice(variableAlignment[pops[[2]]], biAllelicSites.var)
                               p3Slice <- populationSlice(variableAlignment[pops[[3]]], biAllelicSites.var)
                               p4Slice <- populationSlice(variableAlignment[pops[[4]]], biAllelicSites.var)
                               
                               
                               
                               
                               observedAbbaBabaSites <- getAbbaBabaSites(dna)
                               observedNumberOfAbba <- length(observedAbbaBabaSites[[1]])
                               observedNumberOfBaba <- length(observedAbbaBabaSites[[2]])
                               observedD <- calculateD(observedNumberOfAbba, observedNumberOfBaba)
                               
                               
                               
                               results <<- append(results, ABBABABArecord$new(      ))
                             },
                           ## https://github.com/johnomics/Martin_Davey_Jiggins_evaluating_introgression_statistics/blob/master/compare_f_estimators.r
                           
                           
                           
                           )
                         )

getNumberOfAlleles <- function(dna){
  return(colSums(consensusMatrix(dna) != 0))
}

getBiAllelicSites <- function(alleleNum){
  return(which(alleleNum == 2))
}

getVariableSites <- function(alleleNum){
  return(which(alleleNum > 1))
}

subsetSequence <- function(dna, indexes){
  subSeqs <- DNAStringSet(character(length = length(dna)))
  for(i in 1:length(dna)){
    subSeqs[[i]] <- dna[[i]][indexes]
  }
  return(subSeqs)
}

populationSlice <- function(popSeqs, biAllelicSites){
  counts <- consensusMatrix(popSeqs)
  alleles <- colSums(counts != 0)
  S <- sum(alleles > 1)
  return(list(counts = counts, countsBi = counts[,biAllelicSites], alleles = alleles, S = S))
}





getTwoStates <- function(dna){
  twoStates <- which(colSums(consensusMatrix(dna) != 0) == 2)
  twoStatesSeq <- DNAStringSet(character(length = 4))
  for(i in 1:4){
    twoStatesSeq[[i]] <- dna[[i]][twoStates]
  }
  return(list(Sites = twoStates, Seq = twoStatesSeq))
}

getAbbaBabaStates <- function(dna){
  twoStates <- getTwoStates(dna)
  oneNotTwo <- colSums(consensusMatrix(twoStates[["Seq"]][c(1, 2)]) != 0) > 1
  threeNotFour <- colSums(consensusMatrix(twoStates[["Seq"]][c(3, 4)]) != 0) > 1
  abbaBabbaSites <- which(oneNotTwo & threeNotFour)
  abbaBabbaSeq <- DNAStringSet(character(length = 4))
  for(i in 1:4){
    abbaBabbaSeq[[i]] <- twoStatesSeq[[i]][abbaBabbaSites]
  }
  abbaSites <- which(colSums(consensusMatrix(abbaBabbaSeq[c(2,3)]) != 0) == 1)
  babaSites <- which(colSums(consensusMatrix(abbaBabbaSeq[c(1,3)]) != 0) == 1)
  whichSitesAbba <- twoStates[abbaBabbaSites[abbaSites]]
  whichSitesBaba <- twoStates[abbaBabbaSites[babaSites]]
  return(list(AbbaSites = whichSitesAbba, BabaSites = whichSitesBaba))
}





calculateGreenS <- function(abba, babba){
  return(abba - babba)
}

calculateGreenD <- function(abba, baba){
  return(calculateGreenS(abba, baba) / (abba + baba))
}