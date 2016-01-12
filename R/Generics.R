# Manipulating Genetic data.

#' Get the length of the aligned sequences.
#' 
#' @param object The object storing the multiple sequence alignment data from
#' which you wish to get the length.
#' 
#' @return A numeric value representing the length of the aligned sequences.
#' All sequences in the alignment have the same length.
#' @rdname sequenceLength
setGeneric("sequenceLength", function(object){
  standardGeneric("sequenceLength")
})

setGeneric("maskConservedSites", function(object, append = "union"){
  standardGeneric("maskConservedSites")
})

setGeneric("maskNs", function(object, append = "union"){
  standardGeneric("maskNs")
})

setGeneric("maskedSites", function(object){
  standardGeneric("maskedSites")
})

setGeneric("unmaskedSites", function(object){
  standardGeneric("unmaskedSites")
})

setGeneric("nMaskedSites", function(object){
  standardGeneric("nMaskedSites")
})

setGeneric("nUnMaskedSites", function(object){
  standardGeneric("nUnMaskedSites")
})

setGeneric("maskSequences", function(object, seqnames, invert = TRUE, append = "union"){
  standardGeneric("maskSequences")
})

setGeneric("getSeqNames", function(object){
  standardGeneric("getSeqNames")
})

setGeneric("statesPerBase", function(object){
  standardGeneric("statesPerBase")
})

setGeneric("polySites", function(object){
  standardGeneric("polySites")
})

setGeneric("slidingPoly", function(object, windowSize, stepSize){
  standardGeneric("slidingPoly")
})

setGeneric("getWindowsOver", function(object){
  standardGeneric("getWindowsOver")
})