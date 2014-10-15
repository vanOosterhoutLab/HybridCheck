#' @title ComparrisonSettings Reference Class.
#' @description A Reference Class to represent a settings for Triplets.
#'
#' @field Method An integer vector of length 1.
#' @import methods
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
                                         if(length(value != 1)){stop("You must provide one boolean value for the Refine parameter.")}
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
                                     
                                     showSettings =
                                       function(){
                                         message("Settings for Generation of Sequence Scan Combinations.")
                                         return(paste('Triplet Generation Method (Method): ', getMethod(),
                                                '\n\nDistance Threshold for eliminating triplets with\n\ttoo similar sequences (DistanceThreshold): ', getDistanceThreshold(),
                                                '\n\nHow many sequences from the same partition are\n\tallowed in a triplet (PartitionStrictness): ', getPartitionStrictness(),
                                                '\n\nRefine triplet selection (TRUE), or start with\n\ta full set of all possible triplets again (Refine): ', getRefine(), sep=""))
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
  return(class(value) == "numeric" && length(value) == 3)
}