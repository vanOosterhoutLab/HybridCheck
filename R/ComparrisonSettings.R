#' A reference class to represent settings for triplet generation.
#' @name ComparrisonSettings.
#' @field Method An integer vector of length 1.
#' @field DistanceThreshold A numeric vector of length 1.
#' @field PartitionStrictness An integer vector of length 1.
#' @field Refine A logical vector of length 1.
ComparrisonSettings <- setRefClass("ComparrisonSettings",
                                   
                                   fields = list(
                                     SeqNames = "character",
                                     Method = "integer",
                                     DistanceThreshold = "numeric",
                                     PartitionStrictness = "integer",
                                     partiallySignificant = "logical",
                                     TripletCombinations = "list",
                                     AcceptedCombinations = "list"),
                                   
                                   methods = list(
                                     initialize = 
                                       function(dna, ftt){
                                         "Creates the object and sets all parameters to their default."
                                         Method <<- 1L
                                         DistanceThreshold <<- 0.01
                                         PartitionStrictness <<- 2L
                                         partiallySignificant <<- FALSE
                                         SeqNames <<- dna$getSequenceNames()
                                         TripletCombinations <<- combn(dna$getSequenceNames(), 3, simplify=FALSE)
                                         AcceptedCombinations <<- list()
                                         decideAcceptedTriplets(dna, ftt)
                                       },
                                     
                                     changeMethod =
                                       function(value){
                                         value <- unique(value)
                                         if(!is.integer(value)){stop("You must enter integer values between 1 and 4.")}
                                         if(any(value > 4L)){stop("4L is the maximum value allowed.")}
                                         if(any(value < 1L)){stop("1L is the minimum value allowed.")}
                                         if(3L %in% value && 4L %in% value){stop("Choose either method 3L or 4L, not both.")}
                                         Method <<- value
                                       },
                                     
                                     setDistanceThreshold =
                                       function(value){
                                         if(length(value) != 1 || !(value > 0 && value < 1)){stop("You must give a single double value, between 0 and 1, as a distance threshold.")}
                                         DistanceThreshold <<- value
                                       },
                                     
                                     setPartitionStrictness =
                                       function(value){
                                         if(length(value) != 1 || !any(value == c(1L, 2L))){stop("You must enter a single integer value of 1 or 2 as a PartitionStrictness.")}
                                         PartitionStrictness <<- value
                                       },
                                     
                                     setPartiallySignificant =
                                       function(value){
                                         if(length(value) > 1){stop("Only provide one logical value.")}
                                         partiallySignificant <<- value
                                       },
                                     
                                     decideAcceptedTriplets =
                                       function(dna, ftt){
                                         AcceptedCombinations <<- TripletCombinations
                                         rejects <- c()
                                         if(hasMultipleCombinations()){
                                           if(1L %in% Method){
                                             message("Deciding triplets based on results of Four Taxon Tests.")
                                             if(partiallySignificant){
                                               rejects <- c(rejects, fttDescision(ftt, "PART.SIGNIFICANT", TripletCombinations, dna))
                                             } else {
                                               message("Using tests that are globally significant.")
                                               rejects <- c(rejects, fttDescision(ftt, "SIGNIFICANT", TripletCombinations, dna))
                                             }
                                           }
                                           if(2L %in% Method && dna$hasPopulations()){
                                             rejects <- c(rejects, groupDescision(dna$Populations, TripletCombinations, PartitionStrictness)) 
                                           }
                                           if(3L %in% Method || 4L %in% Method){                                                    
                                             rejects <- c(rejects, distanceDescision(dna, Method, DistanceThreshold, TripletCombinations))
                                           }
                                         } else {
                                           warning("There is only one comparrison possible - presumably only 3 sequences are present.")
                                         }
                                         if(!is.null(rejects) && length(rejects) > 0){
                                           AcceptedCombinations <<- AcceptedCombinations[-rejects]
                                         }
                                       },
                                     
                                     hasTripletCombinations =
                                       function(){
                                         return(length(TripletCombinations) > 0)
                                       },
                                     
                                     hasMultipleCombinations =
                                       function(){
                                         return(length(TripletCombinations) > 1)
                                       },
                                     
                                     numberOfTripletCombinations =
                                       function(){
                                         return(length(TripletCombinations))
                                       },
                                     
                                     numberOfAcceptedCombinations =
                                       function(){
                                         return(length(AcceptedCombinations))
                                       },
                                     
                                     textSummary =
                                       function(){
                                         "Creates a character vector of a summary of the comparrison settings."
                                         return(paste('Settings for Sequence Scan Combinations:\n',
                                                      '----------------------------------------\n',
                                                      'Triplet Generation Method (Method): ', paste0(Method, collapse = ", "),
                                                      '\n\nDistance Threshold for excluding triplets with\n\ttoo similar sequences (DistanceThreshold): ', DistanceThreshold,
                                                      '\n\nHow many sequences from the same partition are\n\tallowed in a triplet (PartitionStrictness): ', PartitionStrictness, 
                                                      '\n\nA total of ', numberOfAcceptedCombinations(), ' triplets will be compared.', 
                                                      sep=""))
                                       },
                                     
                                     printCombinations =
                                       function(){
                                         return(paste0(lapply(AcceptedCombinations, function(x) paste0(x, collapse=", ")), collapse="\n"))
                                       },
                                     
                                     htmlCombinations =
                                       function(){
                                         return(paste0(lapply(AcceptedCombinations, function(x) paste0(x, collapse=", ")), collapse="<br>"))
                                       },
                                     
                                     show = function(){
                                       "Prints a summary of the object to console."
                                       cat(textSummary())
                                     },
                                     
                                     setSettings =
                                       function(dna, ftt, ...){
                                         settings <- list(...)
                                         parameters <- names(settings)
                                         for(i in 1:length(settings)){
                                           if(parameters[i] == "Method"){
                                             changeMethod(settings[[i]])
                                           }
                                           if(parameters[i] == "DistanceThreshold"){
                                             setDistanceThreshold(settings[[i]])
                                           }
                                           if(parameters[i] == "PartitionStrictness"){
                                             setPartitionStrictness(settings[[i]])
                                           }
                                           if(parameters[i] == "partiallySignificant"){
                                             setPartiallySignificant(settings[[i]])
                                           }
                                         }
                                         decideAcceptedTriplets(dna, ftt)
                                       }
                                     )
                                   )

#' @name ftt Descision
#' @description Refines the sequence triplets that are to be generated by looking at 
#' which four taxon tests in the FTTmodule were found to be significant.
fttDescision <- function(ftt, significanceStatement, combos, dna){
  significantCombos <- ftt$getFTTs(significanceStatement)
  if(length(significantCombos) > 0){
    popnames <- lapply(significantCombos, function(x) x$getPops())
    matches <- sort(unique(unlist(lapply(popnames, function(x){
      seqnames <- dna$Populations[x]
      which(unlist(lapply(combos, function(y){
        return(sum(any(y %in% seqnames[[1]]), any(y %in% seqnames[[2]]),
                   any(y %in% seqnames[[3]]), any(y %in% seqnames[[4]])) == 3)
        })))
      }))))
    rejects <- which(!1:length(combos) %in% matches)
  } else {
    warning("No four taxon tests were found to be significant, either none were significant or no test has been performed.")
    rejects <- numeric()
  }
  return(rejects)
}

groupDescision <- function(populations, combos, thresh){
  message("Generating triplets to find recombination in sequences, between partitions.")
  matches <- unlist(lapply(populations, function(x){
    which(unlist(lapply(combos, function(y){
      length(which(x %in% y)) > thresh
    })))
  }))
  return(matches)
}

distanceDescision <- function(dna, method, thresh, combos){
  rejects <- c()
  distances <- dist.dna(as.DNAbin(dna$FullSequence), model = "raw")
  seqpairs <- combn(dna$getSequenceNames(), 2, simplify=FALSE)
  if(3L %in% method){
    rejectiondistances <- seqpairs[which(distances < thresh)]
  }
  if(4L %in% method){
    distances_density <- density(distances)
    Lows <- cbind(distances_density$x[which(diff(sign(diff(distances_density$y))) == 2)], distances_density$y[which(diff(sign(diff(distances_density$y))) == 2)])
    Lowest <- Lows[which(Lows[,1] == min(Lows[,1])),]
    rejectiondistances <- seqpairs[which(distances < Lowest[1])]
  }
  for(i in 1:length(combos)){
    for(n in 1:length(rejectiondistances)){
      if(all(rejectiondistances[[n]] %in% combos[[i]])){
        rejects <- c(rejects, i)
        break
      }
    }
  }
  return(rejects)
}
