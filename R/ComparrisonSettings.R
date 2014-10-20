#' A reference class to represent settings for triplet generation.
#' @name ComparrisonSettings.
#' @description
#' @field Method An integer vector of length 1.
#' @field DistanceThreshold A numeric vector of length 1.
#' @field PartitionStrictness An integer vector of length 1.
#' @field Refine A logical vector of length 1.
ComparrisonSettings <- setRefClass("ComparrisonSettings",
                                   
                                   fields = list(
                                     Method = "integer",
                                     DistanceThreshold = "numeric",
                                     PartitionStrictness = "integer",
                                     Refine = "logical",
                                     TripletCombinations = "list"),
                                   
                                   methods = list(
                                     initialize = 
                                       function(){
                                         "Creates the object and sets all parameters to their default."
                                         Method <<- 1L
                                         DistanceThreshold <<- 0.01
                                         PartitionStrictness <<- 2L
                                         Refine <<- FALSE
                                         TripletCombinations <<- list()
                                       },
                                     
                                     getMethod =
                                       function(){
                                         return(Method)
                                       },
                                     
                                     changeMethod =
                                       function(value){
                                         if(length(value) != 1 || !any(value == c(1L, 2L, 3L, 4L))){stop("You must enter a single integer value between 1 and 4.")}
                                         Method <<- value
                                       },
                                     
                                     getDistanceThreshold =
                                       function(){
                                         return(DistanceThreshold)
                                       },
                                     
                                     setDistanceThreshold =
                                       function(value){
                                         if(length(value) != 1 || !(value > 0 && value < 1)){stop("You must give a single double value, between 0 and 1, as a distance threshold.")}
                                         DistanceThreshold <<- value
                                       },
                                     
                                     getPartitionStrictness =
                                       function(){
                                         return(PartitionStrictness)
                                       },
                                     
                                     setPartitionStrictness =
                                       function(value){
                                         if(length(value) != 1 || !any(value == c(1L, 2L))){stop("You must enter a single integer value of 1 or 2 as a PartitionStrictness.")}
                                         PartitionStrictness <<- value
                                       },
                                     
                                     getRefine =
                                       function(){
                                         return(Refine)
                                       },
                                     
                                     setRefine =
                                       function(value){
                                         if(length(value) > 1){stop("You must provide one boolean value for the Refine parameter.")}
                                         Refine <<- value
                                       },
                                     
                                     getTripletCombinations =
                                       function(){
                                         return(TripletCombinations)
                                       },
                                     
                                     setTripletCombinations =
                                       function(value){
                                         if(length(value) < 1 || !is.integertriplet(value)){stop("Triplet Combinations must be a list of numeric vectors, each of length 3.")}
                                         TripletCombinations <<- value
                                       },
                                     
                                     hasTripletCombinations =
                                       function(){
                                         return(length(TripletCombinations) > 0)
                                       },
                                     
                                     hasMultipleCombinations =
                                       function(){
                                         return(length(TripletCombinations) > 1)
                                       },
                                     
                                     numberOfTripletCombos =
                                       function(){
                                         return(length(TripletCombinations))
                                       },
                                     
                                     eliminateTripletCombinations =
                                       function(rejects){
                                         TripletCombinations <<- TripletCombinations[-unlist(rejects)]
                                       },
                                     
                                     textSummary =
                                       function(){
                                         "Creates a character vector of a summary of the comparrison settings."
                                         return(paste('Settings for Sequence Scan Combinations:\n',
                                                      '----------------------------------------\n',
                                                      'Triplet Generation Method (Method): ', getMethod(),
                                                '\n\nDistance Threshold for eliminating triplets with\n\ttoo similar sequences (DistanceThreshold): ', getDistanceThreshold(),
                                                '\n\nHow many sequences from the same partition are\n\tallowed in a triplet (PartitionStrictness): ', getPartitionStrictness(),
                                                '\n\nRefine triplet selection (TRUE), or start with\n\ta full set of all possible triplets again (Refine): ', getRefine(), sep=""))
                                       },
                                     
                                     show = function(){
                                       "Prints a summary of the object to console."
                                       cat(textSummary())
                                     },
                                     
                                     setSettings =
                                       function(...){
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
                                         }
                                       }
                                     
                                     )
                                   )

is.integertriplet <- function(value){
  return(all(unlist(lapply(value, function(x){is.numeric(x) && length(x) == 3}))))
}

generateTriplets2 <- function(settings, ranges){
  message("Generating triplets to find recombination in sequences, between partitions.")
  return(unlist(lapply(ranges, function(x){
    which(unlist(lapply(settings$getTripletCombinations(), function(y){
      length(which(x %in% y)) > settings$getPartitionStrictness()
    })))
  })))
}

generateTripletsDist <- function(dna, settings){
  rejects <- c()
  distances <- dist.dna(as.DNAbin(dna$FullSequence), model = "raw")
  seqpairs <- combn(1:dna$numberOfSequences(), 2, simplify=FALSE)
  if(settings$getMethod() == 3){
    rejectiondistances <- seqpairs[which(distances < settings$getDistanceThreshold())]
  } else {
    rejectiondistances <- generateTriplets4(distances)
  }
  for(i in 1:settings$numberOfTripletCombos()){
    for(n in 1:length(rejectiondistances)){
      if(all(rejectiondistances[[n]] %in% settings$TripletCombinations[[i]])){
        rejects <- append(rejects, i)
        break
      }
    }
  }
  return(rejects)
}

generateTriplets4 <- function(distances){
  distances_density <- density(distances)
  Lows <- cbind(distances_density$x[which(diff(sign(diff(distances_density$y))) == 2)], distances_density$y[which(diff(sign(diff(distances_density$y))) == 2)])
  Lowest <- Lows[which(Lows[,1] == min(Lows[,1])),]
  return(seqpairs[which(distances < Lowest[1])])
}