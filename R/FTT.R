#' A Reference Class for performing and storing results of four taxon tests for introgression.
#' @name FTTester
#' @field global Logical, whether a global statistic will be calculated for the four taxon tests.
#' @field results A list of reference objects defining the result of a given four taxon test.
FTTester <- setRefClass("FTTester",
                         
                         fields = list(blockSize = "integer",
                                       taxaCombos = "list",
                                       results = "list"),
                         
                         methods = list(
                           initialize =
                             function(dna){
                               "Initialization method creates object with its default values."
                               blockSize <<- 50000L
                               results <<- list()
                             },
                           
                           setTaxaCombos =
                             function(taxas, dna){
                               if(length(taxas) > 0){
                                 for(i in taxas){
                                   if(length(unique(i)) != 4){stop("Each taxon combination must provide 4 unique populations: a P1, a P2, a P3 and an A.")}
                                   if(!is.null(names(i))){
                                     inputCheck1 <- all(names(i) %in% c("P1", "P2", "P3", "A"))
                                     if(!inputCheck1){stop("The only names allowed for specifying taxa combos are 'P1', 'P2', 'P3', and 'A'")}
                                   } else {
                                     warning(paste0("No names were provided for the population combination: ", paste0(i, collapse=", "), ".\n",
                                     "Assuming that P1 is ", i[1], ", that P2 is ", i[2], ", that P3 is ", i[3], " and that A is ", i[4]))
                                     names(i) <- c("P1", "P2", "P3", "A")
                                   }
                                   if(!all(i %in% dna$namesOfPopulations())){stop("You have listed a population name or number in a taxa combo which does not exist.")}
                                   taxaCombos <<- c(taxaCombos, list(i))
                                 }
                               }
                             },
                           
                           hasTaxaCombos =
                             function(){
                               return(length(taxaCombos) > 0)
                             },
                           
                           setBlockSize =
                             function(value){
                               if(!is.integer(value) || length(value != 1)){stop("Provide only one, integer value.")}
                               blockSize <<- value
                             },
                          
                           setSettings =
                             function(...){
                               settings <- list(...)
                               parameters <- names(settings)
                               for(i in 1:length(settings)){
                                 if(parameters[i] == "blockSize"){
                                   setBlockSize(settings[[i]])
                                 }
                                 if(parameters[i] == "taxaCombos"){
                                   setTaxaCombos(settings[[i]])
                                 }
                               }
                             },
                           
                           generateFTTs =
                             function(){
                               message("Initializing new FTtest data.")
                               results <<- lapply(taxaCombos, function(x) FTTrecord$new(x$P1, x$P2, x$P3, x$A))
                             },
                           
                           runFTTs =
                             function(){
                               
                             },
                           
                           calculateDandFd =
                             function(aln, pops){
                               # State counts at each site.
                               counts.all <- consensusMatrix(aln)
                               # Calculate the number of alleles at each site in the alignment.
                               numAlleles.all <- colSums(counts.all != 0)
                               # Find which sites are bi-allelic.
                               biSites.all <- which(numAlleles.all == 2)
                               # Find which ones are variable.
                               varSites.all <- which(numAlleles.all > 1)
                               # Make a version of the sequence alignment which only includes variable sites.
                               aln.var <- subsetSequence(aln, varSites.all)
                               # Get the number of alleles at those variable sites.
                               numAlleles.var <- numAlleles.all[varSites.all]
                               S.all <- length(numAlleles.var)
                               # Find which of the variable sites are bi-allelic.
                               biSites.var <- which(numAlleles.var == 2)
                               # Population slices - refer to function's description.
                               p1Slice <- populationSlice(aln.var[pops[[1]]], biSites.var)
                               p2Slice <- populationSlice(aln.var[pops[[2]]], biSites.var)
                               p3Slice <- populationSlice(aln.var[pops[[3]]], biSites.var)
                               p4Slice <- populationSlice(aln.var[pops[[4]]], biSites.var)
                               statsResults <- calculateStats(numberOfAlleles.all, biAllelicSites.all,
                                                              p1Slice, p2Slice, p3Slice, p4Slice) 
                             }
                           )
                         )


FTTrecord <- setRefClass("FTTrecord",
                         
                         fields = list(
                           P1 = "character",
                           P2 = "character",
                           P3 = "character",
                           A = "character",
                           blockSize = "integer",
                           recordFile = "character",
                           table = "data.frame"
                           ),
                         
                         methods = list(
                           initialize =
                             function(p1, p2, p3, a, bs){
                               P1 <<- p1
                               P2 <<- p2
                               P3 <<- p3
                               A <<- a
                               blockSize <<- bs
                               
                             }
                           )
                         )


#' @name Subset a DNAStringSet.
#' @description For extracting only certain cites of an MSA represented as a DNA
subsetSequence <- function(dna, indexes){
  subSeqs <- DNAStringSet(character(length = length(dna)))
  for(i in 1:length(dna)){
    subSeqs[[i]] <- dna[[i]][indexes]
  }
  return(subSeqs)
}

#' @name populationSlice
#' @description Calculate for a subset of an alignment of variable sites:
#' The counts, the number of alleles, and the counts of the (global) bi-allelic,
#' sites in the population. 
populationSlice <- function(popSeqs, biSites){
  counts <- consensusMatrix(popSeqs)
  alleles <- colSums(counts != 0)
  S <- sum(alleles > 1)
  return(list(counts = counts, countsBi = counts[, biSites], alleles = alleles, S = S))
}

calculateStats <- function(counts.all, biSites.all, slice1, slice2, slice3, slice4){
  ABBA <- 0
  BABA <- 0
  maxABBA_D <- 0
  maxBABA_D <- 0
  for(i in 1:length(biSites.all)){
    i.bi.counts <- counts.all[, biSites.all[i]]
    alleles <- names(which(i.bi.counts > 0))
    anc <- names(which(slice4$countsBi[, i] != 0))
    if(length(anc) != 1){
      anc <- names(which(i.bi.counts == max(i.bi.counts)))[1]
    }
    derived <- alleles[which(alleles != anc)]
    P1df <- slice1$countsBi[derived, i] / sum(slice1$countsBi[, i])
    P2df <- slice2$countsBi[derived, i] / sum(slice2$countsBi[, i])
    P3df <- slice3$countsBi[derived, i] / sum(slice3$countsBi[, i])
    P4df <- slice4$countsBi[derived, i] / sum(slice4$countsBi[, i])
    ABBA <- sum(ABBA, (1 - P1df) * P2df * P3df * (1 - P4df), na.rm = TRUE)
    BABA <- sum(BABA, P1df * (1 - P2df) * P3df * (1 - P4df), na.rm = TRUE)
    if (is.na(P3df) == FALSE & is.na(P2df) == FALSE & P3df >= P2df){
      maxABBA_D <- sum(maxABBA_D, (1 - P1df) * P3df * P3df * (1 - P4df), na.rm = TRUE)
      maxBABA_D <- sum(maxBABA_D, P1df * (1 - P3df) * P3df * (1 - P4df), na.rm = TRUE)
    } else {
      maxABBA_D <- sum(maxABBA_D, (1 - P1df) * P2df * P2df * (1 - P4df), na.rm = TRUE)
      maxBABA_D <- sum(maxBABA_D, P1df * (1 - P2df) * P2df * (1 - P4df), na.rm = TRUE)
    }
  }
  out <- c(ABBA=ABBA, BABA=BABA,
           D = (ABBA - BABA) / (ABBA + BABA),
           maxABBA_D = maxABBA_D, maxBABA_D = maxBABA_D,
           Fd = (ABBA - BABA) / (maxABBA_D - maxBABA_D))
  return(out)
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