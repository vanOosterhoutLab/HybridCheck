#' @include Generic_Methods.R
NULL

setMethod("sequenceLength",
          representation("DNAStringSet"),
          function(object) {
            return(width(object))
          })


setMethod("statesPerBase",
          signature("DNAStringSet"),
          function(object) {
            return(colSums(consensusMatrix(object) != 0))
          })

setMethod("polymorphicSites",
          signature("DNAStringSet"),
          function(object) {
            return(statesPerBase(object) > 1)
          })

setMethod("sitesWithUnknown",
          signature("DNAStringSet"),
          function(object){
            consensusMatrix(object)[15, ] != 0
          })


