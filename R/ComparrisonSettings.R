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
                                     Groups = "list",
                                     TripletCombinations = "list",
                                     AcceptedCombinations = "list"),
                                   
                                   methods = list(
                                     initialize = 
                                       function(dna){
                                         "Creates the object and sets all parameters to their default."
                                         Method <<- 1L
                                         DistanceThreshold <<- 0.01
                                         PartitionStrictness <<- 2L
                                         SeqNames <<- dna$getSequenceNames()
                                         TripletCombinations <<- combn(dna$getSequenceNames(), 3, simplify=FALSE)
                                         AcceptedCombinations <<- list()
                                         decideAcceptedTriplets()
                                       },
                                     
                                     changeMethod =
                                       function(value){
                                         value <- unique(value)
                                         if(!is.integer(value)){stop("You must enter integer values between 1 and 3.")}
                                         if(any(value > 3L)){stop("3L is the maximum value allowed.")}
                                         if(any(value < 1L)){stop("1L is the minimum value allowed.")}
                                         if(2L %in% value && 3L %in% value){stop("Select either method 2 or 3, not both.")}
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
                                     
                                     setGroups =
                                       function(groups){
                                         if(length(groups) > 0){
                                           if(any(!unlist(lapply(groups, function(x) is.integer(x) || is.character(x))))){stop("Need to provide a list of groups of sequence names or integers representing sequence numbers.")}
                                           groups <- lapply(groups, function(x){
                                             if(is.numeric(x)){
                                               return(SeqNames[x]) 
                                             } else {
                                               return(x)
                                             }
                                           })
                                           if(any(table(unlist(groups)) > 1)){stop("Entered a sequence name or number in more than one group.")}
                                           if(any(!unlist(lapply(groups, function(x) all(x %in% SeqNames))))){stop("Some sequences specified in the groups are not in the sequence data.")}
                                         }
                                         Groups <<- groups
                                         },
                                     
                                     decideAcceptedTriplets =
                                       function(dna){
                                         AcceptedCombinations <<- TripletCombinations
                                         rejects <- c()
                                         if(hasMultipleCombinations()){
                                           if(1 %in% Method){
                                             rejects <- c(rejects, groupDescision(Groups, TripletCombinations, PartitionStrictness)) 
                                           }
                                           if(2 %in% Method || 3 %in% Method){                                                    
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
                                                      'Triplet Generation Method (Method): ', Method,
                                                      '\n\nSequences are organized according to the following groups: \n',
                                                      paste(Groups, collapse = ',\n'),
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
                                       function(dna, ...){
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
                                           if(parameters[i] == "Refine"){
                                             setRefine(settings[[i]])
                                           }
                                           if(parameters[i] == "Groups"){
                                             setGroups(settings[[i]])
                                           }
                                         }
                                         decideAcceptedTriplets(dna)
                                       }
                                     )
                                   )

groupDescision <- function(groups, combos, thresh){
  message("Generating triplets to find recombination in sequences, between partitions.")
  matches <- unlist(lapply(groups, function(x){
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
  if(2L %in% method){
    rejectiondistances <- seqpairs[which(distances < thresh)]
  }
  if(3L %in% method){
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
