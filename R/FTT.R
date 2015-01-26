#' A Reference Class for performing and storing results of four taxon tests for introgression.
#' @name FTTester
#' @field global Logical, whether a global statistic will be calculated for the four taxon tests.
#' @field results A list of reference objects defining the result of a given four taxon test.
FTTester <- setRefClass("FTTester",
                         
                         fields = list(numBlocks = "integer",
                                       taxaCombos = "list",
                                       results = "list"),
                         
                         methods = list(
                           initialize =
                             function(dna){
                               "Initialization method creates object with its default values."
                               numBlocks <<- 50000L
                               results <<- list()
                             },
                           
                           manualTaxaCombos =
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
                                   if(!all(i %in% dna$namesOfPopulations())){stop("You have listed a population name in a taxa combo which does not exist.")}
                                   taxaCombos <<- c(taxaCombos, list(i))
                                 }
                               }
                             },
                           
                           autoTaxaCombos =
                             function(dna){
                               if(dna$numberOfPopulations() < 4){
                                 stop("Less than 4 populations have been defined - can't form a quartet for an ABBA-BABA test.")
                               }
                               allCombs <- combn(dna$namesOfPopulations(), 4, simplify = FALSE)
                               allDists <- as.matrix(stringDist(dna$FullSequence, method = "hamming")) / dna$getFullLength()
                               generatedCombs <- lapply(allCombs, function(x){
                                 out <- list(P1 = NULL, P2 = NULL, P3 = NULL, A = NULL)
                                 otus <- x
                                 seqsInOtus <- dna$Populations[otus]
                                 otuPairs <- combn(otus, 2, simplify = FALSE)
                                 distances <- compDist(otuPairs, seqsInOtus, allDists)
                                 minOTUs <- distances[which(distances$dist == min(distances$dist)), c("OTU1", "OTU2")]
                                 out$P1 <- minOTUs$OTU1
                                 out$P2 <- minOTUs$OTU2
                                 remainingOtus <- otus[which(otus != out$P1 & otus != out$P2)]
                                 seqsInOtus2 <- dna$Populations[remainingOtus]
                                 P1P2otu <- paste0("(", out$P1, ", ", out$P2, ")")
                                 remainingOtus <- c(remainingOtus, P1P2otu)
                                 seqsInOtus2[[P1P2otu]] <- unlist(seqsInOtus[c(out$P1, out$P2)])
                                 remainingOtuPairs <- combn(remainingOtus, 2, simplify = FALSE)
                                 remainingOtuPairs <- remainingOtuPairs[unlist(lapply(remainingOtuPairs, function(x) P1P2otu %in% x))]
                                 distances_2 <- compDist(remainingOtuPairs, seqsInOtus2, allDists)
                                 minOTUs_2 <- distances_2[which(distances_2$dist == min(distances_2$dist)), c("OTU1", "OTU2")]
                                 out$P3 <- minOTUs_2[1, which(minOTUs_2 != P1P2otu)]
                                 out$A <- otus[which(!(otus %in% unlist(out)))]
                                 checks <- c(
                                   # Check P1 & P2 are closer than P1 and P3.
                                   "d(P1,P2) < d(P1,P3)" = distances[which(apply(distances, 1, function(x){(out$P1 %in% x) & (out$P2 %in% x)})), 3] < 
                                     distances[which(apply(distances, 1, function(x) (out$P1 %in% x) & (out$P3 %in% x))), 3],
                                   # Check P1 & P2 are closer than P1 and A.
                                   "d(P1,P2) < d(P1,PA)" = distances[which(apply(distances, 1, function(x) (out$P1 %in% x) & (out$P2 %in% x))), 3] < 
                                     distances[which(apply(distances, 1, function(x) (out$P1 %in% x) & (out$A %in% x))), 3],
                                   # Check P1 & P2 are closer than P2 and P3.
                                   "d(P1,P2) < d(P2,P3)" = distances[which(apply(distances, 1, function(x){(out$P1 %in% x) & (out$P2 %in% x)})), 3] < 
                                     distances[which(apply(distances, 1, function(x) (out$P2 %in% x) & (out$P3 %in% x))), 3],
                                   # Check P1 & P2 are closer than P2 and A.
                                   "d(P1,P2) < d(P2,A)" = distances[which(apply(distances, 1, function(x) (out$P1 %in% x) & (out$P2 %in% x))), 3] < 
                                     distances[which(apply(distances, 1, function(x) (out$P2 %in% x) & (out$A %in% x))), 3]
                                 )
                                 if(any(!checks)){
                                   warning(paste0(paste0("Automatically allocating P1, P2, P3 and A designations for population quartet ", paste0(otus, collapse=", "),":\n"),
                                                  paste0("Sanity check: ", names(checks)[which(!checks)], " failed.", collapse = "\n")))
                                 }
                                 return(out)
                               })
                               manualTaxaCombos(generatedCombs, dna)
                             },
                           
                           hasTaxaCombos =
                             function(){
                               return(length(taxaCombos) > 0)
                             },
                           
                           setNumBlocks =
                             function(value){
                               if(!is.integer(value) || length(value != 1)){stop("Provide only one, integer value.")}
                               numBlocks <<- value
                             },
                          
                           setSettings =
                             function(...){
                               settings <- list(...)
                               parameters <- names(settings)
                               for(i in 1:length(settings)){
                                 if(parameters[i] == "numBlocks"){
                                   setNumBlocks(settings[[i]])
                                 }
                                 if(parameters[i] == "taxaCombos"){
                                   setTaxaCombos(settings[[i]])
                                 }
                               }
                             },
                           
                           generateFTTs =
                             function(hybridsDir){
                               message("Initializing new FTtest data.")
                               results <<- lapply(taxaCombos, function(x) FTTrecord$new(x[["P1"]], x[["P2"]], x[["P3"]], x[["A"]], hybridsDir))
                             },
                           
                           getAllNames = 
                             function(){
                               return(lapply(results, function(x) x$getPops()))
                             },
                           
                           printAllNames =
                             function(){
                               quadNames <- lapply(getAllNames(), function(x) paste(x, collapse = ", "))
                               collectedNames <- lapply(1:length(quadNames), function(i) paste(i, quadNames[[i]], sep=": "))
                               return(collectedNames)
                             },
                           
                           getFTTs =
                             function(selections){
                               "Returns a list of references to FTTrecords objects according to user selection."
                               if(!is.list(selections)){
                                 selections <- list(selections)
                               }
                               selections <- unique(selections)
                               if(!is.null(selections) && length(selections) > 0){
                                 ind <- numeric()
                                 if(length(results) != 0){
                                   if("ALL" %in% selections){
                                     ind <- 1:length(results)
                                   } else {
                                     if("NOT.TESTED" %in% selections){
                                       ind <- c(ind, which(unlist(lapply(results, function(x) x$noTestPerformed()))))
                                       selections <- selections[which(selections != "NOT.TESTED")]
                                     }
                                     if("TESTED" %in% selections){
                                       ind <- c(ind, which(unlist(lapply(results, function(x) !x$noTestPerformed()))))
                                       selections <- selections[which(selections != "TESTED")]
                                     }
                                     if("SIGNIFICANT" %in% selections){
                                       ind <- c(ind, which(unlist(lapply(results, function(x) x$globallySignificant()))))
                                       selections <- selections[which(selections != "SIGNIFICANT")]
                                     }
                                     if("PART.SIGNIFICANT" %in% selections){
                                       ind <- c(ind, which(unlist(lapply(results, function(x) x$globallySignificant()))))
                                       selections <- selections[which(selections != "PART.SIGNIFICANT")]
                                     }
                                     if(any(unlist(lapply(selections, length)) != 4)){stop("Selections must provide a vector of 4 sequence names.")}
                                     if(any(unlist(lapply(selections, function(x) !is.character(x))))){stop("Selections must be of class character.")}
                                     allNames <- do.call(rbind, getAllNames())
                                     ind <- c(ind, unlist(lapply(selections, function(x){
                                       which(allNames[,1] %in% x & allNames[,2] %in% x & allNames[,3] %in% x & allNames[,4] %in% x)
                                     })))
                                   }
                                 }
                                 ind <- unique(ind)
                                 return(results[ind])
                               } else {
                                 stop("No selection of FTT was provided.")
                               }
                             },
                           
                           runFTTests =
                             function(selections, dna, numBlocks, blocksLen){
                               fttsToTest <- getFTTs(selections)
                               for(ftt in fttsToTest){
                                 fourTaxonTest(dna, ftt, numBlocks, blocksLen)
                               }
                             },
                           
                           getResults =
                             function(selections){
                               fttsToCollect <- getFTTs(selections)
                               collectedTables <- lapply(fttsToCollect, function(ftt) ftt$getTable())
                               return(do.call(rbind, collectedTables))
                             }
                           )
                         )

#' @name compDist
#' @description Compute hamming distances between a set of composite OTUs. 
compDist <- function(popPairs, seqInPop, distMat){
  distances <- lapply(popPairs, function(y){
    grid <- expand.grid(seqInPop[y])
    ds <- apply(grid, 1, function(z){distMat[z[1], z[2]]})
    return(sum(ds) / nrow(grid))
  })
  out <- data.frame(do.call(rbind, popPairs))
  colnames(out) <- c("OTU1", "OTU2")
  out$OTU1 <- as.character(out$OTU1)
  out$OTU2 <- as.character(out$OTU2)
  out$dist <- unlist(distances)
  return(out)
}

#' @name Subset a DNAStringSet.
#' @description For extracting only certain cites of an MSA represented as a DNA
subsetSequence <- function(dna, indexes){
  subSeqs <- DNAStringSet(character(length = length(dna)))
  for(i in 1:length(dna)){
    subSeqs[[i]] <- dna[[i]][indexes]
  }
  names(subSeqs) <- names(dna) 
  return(subSeqs)
}

calculateStats <- function(counts.all, biSites.all, slice1, slice2, slice3, slice4){
  # Allocate space for the proportion of derived allele in each population.
#   P1df <- vector(mode = "numeric", length = length(biSites.all))
#   P2df <- vector(mode = "numeric", length = length(biSites.all))
#   P3df <- vector(mode = "numeric", length = length(biSites.all))
#   P4df <- vector(mode = "numeric", length = length(biSites.all))
  # Allocate space for Cabba and Cbaba.
  ABBA <- vector(mode = "numeric", length = length(biSites.all))
  BABA <- vector(mode = "numeric", length = length(biSites.all))
  # Allocate space for 
  maxABBA_D <- vector(mode = "numeric", length = length(biSites.all))
  maxBABA_D <- vector(mode = "numeric", length = length(biSites.all))
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
    ABBA[i] <- (1 - P1df) * P2df * P3df * (1 - P4df)
    BABA[i] <- P1df * (1 - P2df) * P3df * (1 - P4df)
    if(!is.na(P3df) & !is.na(P2df) & P3df >= P2df){
      maxABBA_D[i] <- (1 - P1df) * P3df * P3df * (1 - P4df)
      maxBABA_D[i] <- P1df * (1 - P3df) * P3df * (1 - P4df)
    } else {
      maxABBA_D[i] <- (1 - P1df) * P2df * P2df * (1 - P4df)
      maxBABA_D[i] <- P1df * (1 - P2df) * P2df * (1 - P4df)
    }
  }
  numABBA <- length(which(ABBA > BABA))
  numBABA <- length(which(ABBA < BABA))
  if(numABBA < numBABA){
    binomialP <- pbinom(numABBA, (numABBA + numBABA), .5, lower.tail = TRUE)
  } else {
    binomialP <- pbinom(numBABA, (numABBA + numBABA), .5, lower.tail = TRUE)
  }
  out <- data.frame(numABBA = numABBA, numBABA = numBABA, ABBA = sum(ABBA), BABA = sum(BABA), D_P_binom = binomialP,
           maxABBA_D = sum(maxABBA_D), maxBABA_D = sum(maxBABA_D))
  return(out)
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

calculateDandFd <- function(aln, pops){
  # State counts at each site.
  counts.all <- consensusMatrix(aln)
  # Calculate the number of alleles at each site in the alignment.
  num.alleles.all <- colSums(counts.all != 0)
  # Find which sites are bi-allelic.
  bi.all <- which(num.alleles.all == 2)
  # Find which ones are variable.
  var.sites.all <- which(num.alleles.all > 1)
  # Make a version of the sequence alignment which only includes variable sites.
  aln.var <- subsetSequence(aln, var.sites.all)
  # Get the number of alleles at those variable sites.
  alleles.var <- num.alleles.all[var.sites.all]
  s.all <- length(alleles.var)
  # Find which of the variable sites are bi-allelic.
  bi.var <- which(alleles.var == 2)
  # Population slices - refer to function's description.
  p1Slice <- populationSlice(aln.var[pops[[1]]], bi.var)
  p2Slice <- populationSlice(aln.var[pops[[2]]], bi.var)
  p3Slice <- populationSlice(aln.var[pops[[3]]], bi.var)
  p4Slice <- populationSlice(aln.var[pops[[4]]], bi.var)
  return(data.frame(S = s.all, P1_S = p1Slice$S, P2_S = p2Slice$S,
             P3_S = p3Slice$S, P4_S = p4Slice$S, 
             calculateStats(counts.all, bi.all,
                            p1Slice, p2Slice, p3Slice, p4Slice)))
}

#' @title fourTaxonTest
#' @name fourTaxonTest
#' @description Computes the Four Taxon Test for a given set of P1, P2, P3, and P4.
#' Patterson's D is calculated, as is Fd from the paper Martin et al. (2014).
fourTaxonTest <- function(dna, fttRecord, numBlocks, lengthOfBlocks){
  # The jack-knife implementation and computation of Z score has been based on that
  # used in the software ANGSD.
  if(is.null(numBlocks) && is.null(lengthOfBlocks)){
    stop("Invalid input - no number of blocks or size of blocks provided.")
  }
  if(!is.null(numBlocks) && !is.null(lengthOfBlocks)){
    stop("Provide the number of blocks to divide sequence alignment into, OR the proposed size of the blocks, not both.")
  }
  dnaLen <- dna$getFullLength()
  if(!is.null(lengthOfBlocks) && is.null(numBlocks)){
    fttRecord$numBlocks <- as.integer(floor(dnaLen / lengthOfBlocks))
  } else {
    fttRecord$numBlocks <- numBlocks
  }
  fttRecord$blockLength <- as.integer(floor(dnaLen / fttRecord$numBlocks))
  blockStart <- seq(from = 1, to = dnaLen, by = fttRecord$blockLength)
  blockEnd <- seq(from = fttRecord$blockLength, to = dnaLen, by = fttRecord$blockLength)
  if(length(blockStart) > length(blockEnd)){
    blockStart <- blockStart[1:length(blockStart) - abs(length(blockStart) - length(blockEnd))]
  }
  results <- data.frame(BlockStart = blockStart, BlockEnd = blockEnd)
  blocks <- apply(results, 1, function(x){subsetSequence(dna$FullSequence, x[1]:x[2])})
  blocksStats <- do.call(rbind, lapply(blocks, function(x){calculateDandFd(x, dna$Populations[c(fttRecord$P1, fttRecord$P2, fttRecord$P3, fttRecord$A)])}))
  blocksStats$S1234 <- blocksStats$ABBA - blocksStats$BABA
  blocksStats$abbaBabaSum <- blocksStats$ABBA + blocksStats$BABA
  blocksStats$D <- blocksStats$S1234 / blocksStats$abbaBabaSum
  
  blocksStats$S1DD4 <- blocksStats$maxABBA_D - blocksStats$maxBABA_D
  blocksStats$Fd <- blocksStats$S1234 / blocksStats$S1DD4
  blocksStats$fdD0 <- as.numeric(blocksStats$Fd)
  blocksStats$fdD0[c(which(is.na(blocksStats$D)), which(blocksStats$D < 0))] <- NA
  
  # Jack-knifing based on ANGSD implementation.
  # Number of blocks.
  nBlocks <- nrow(blocksStats)
  # Calculates D. 
  dCalc <- function(x){sum(x$S)/sum(x$abbaBabaSum)}
  # Global Estimate.
  globalEstimate <- dCalc(blocksStats[, c("S", "abbaBabaSum")])
  # Pseudoestimates.
  blocksStats$pseudoEstimates <- rep(0, nBlocks)
  blocksStats$blockFraction <- blocksStats$abbaBabaSum / sum(blocksStats$abbaBabaSum)
  for(i in 1:nBlocks){
    blocksStats$pseudoEstimates[i] <- dCalc(blocksStats[-i, c("S1234", "abbaBabaSum")])
  }
  blocksStats$invBlockFraction <- 1 - blocksStats$blockFraction
  blocksStats$scaledPseudoEstimate <- blocksStats$invBlockFraction * blocksStats$pseudoEstimates
  blocksStats$scaledGlobalEstimate <- blocksStats$blockFraction * globalEstimate
  prodGlobN <- nBlocks * globalEstimate
  scaledPseudoEstimateSum <- sum(blocksStats$scaledPseudoEstimate)
  fracReciprocal <- 1 / blocksStats$blockFraction - 1
  fttRecord$D_jEstimate <- prodGlobN - scaledPseudoEstimateSum
  fttRecord$D_jVariance <- 1/nBlocks * 
    sum(1 / fracReciprocal * ( ((1 / blocksStats$blockFraction) * globalEstimate) - 
                                 (fracReciprocal * blocksStats$pseudoEstimates) - (nBlocks * globalEstimate) + scaledPseudoEstimateSum)^2)
  fttRecord$D_jSD <- sqrt(fttRecord$D_jVariance)
  fttRecord$D_jZ <- fttRecord$D_jEstimate / fttRecord$D_jSD
  
  
  
  fttRecord$table <- cbind(results, blocksStats)
  fttRecord$D_globalX2 <- -2 * sum(log(fttRecord$table$D_P_binom))
  fttRecord$fisherBinomialSum <- pchisq(fttRecord$D_globalX2,
                              df = 2 * length(fttRecord$table$D_P_binom), 
                              lower.tail = FALSE)
}

FTTrecord <- setRefClass("FTTrecord",
                         
                         fields = list(
                           P1 = "character",
                           P2 = "character",
                           P3 = "character",
                           A = "character",
                           numBlocks = "integer",
                           blockLength = "integer",
                           D_globalX2 = "numeric",
                           globalABBAcount = "numeric",
                           globalBABAcount = "numeric",
                           globalBinomialP = "numeric",
                           fisherBinomialSum = "numeric",
                           D_jEstimate = "numeric",
                           D_jVariance = "numeric",
                           D_jSD = "numeric",
                           D_jZ = "numeric",
                           
                           tableFile = "character",
                           table = function(value){
                             if(missing(value)){
                               read.table(tableFile)
                             } else {
                               write.table(value, file = tableFile)
                             }
                           }
                           ),
                         
                         methods = list(
                           initialize =
                             function(p1, p2, p3, a, hybridsDir){
                               P1 <<- p1
                               P2 <<- p2
                               P3 <<- p3
                               A <<- a
                               blockLength <<- 0L
                               numBlocks <<- 0L
                               D_globalX2 <<- 0
                               globalBinomialP <<- 0
                               tableFile <<- tempfile(pattern = "FTTtable", tmpdir = hybridsDir)
                               blankTable()
                             },
                           
                           noTestPerformed =
                             function(){
                               return(all(is.na(table)))
                             },
                           
                           globallySignificant =
                             function(){
                               return((!noTestPerformed()) && (globalBinomialP < 0.05))
                             },
                           
                           blankTable =
                             function(){
                               "Method clears the table."
                               table <<- data.frame(BlockStart = NA, BlockEnd = NA,
                                                    ABBA = NA, BABA = NA, maxABBA_D = NA,
                                                    maxBABA_D = NA, D = NA, Fd = NA)
                             },
                           
                           getPops =
                             function(){
                               return(c(P1 = P1, P2 = P2, P3 = P3, A = A))
                             },
                           
                           getTable =
                             function(){
                               return(cbind(P1, P2, P3, A, table, globalX2, globalP))
                             }
                           )
                         )