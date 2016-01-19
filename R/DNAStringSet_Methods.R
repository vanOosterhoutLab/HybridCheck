#' @include Generic_Methods.R
NULL

#' @rdname sequenceLength
setMethod("sequenceLength",
          representation("DNAStringSet"),
          function(object) {
            return(width(object))
          })

#' @rdname statesPerBase
setMethod("statesPerBase",
          signature("DNAStringSet"),
          function(object) {
            return(colSums(consensusMatrix(object) != 0))
          })

#' @rdname polymorphicSites
setMethod("polymorphicSites",
          signature("DNAStringSet"),
          function(object) {
            return(statesPerBase(object) > 1)
          })

#' @rdname sitesHave
setMethod("sitesHave",
          signature("DNAStringSet", "character"),
          function(object, letter){
            mat <- consensusMatrix(object)
            row <- which(rownames(mat) == letter)
            if(length(row) != 1){
              stop("Error subsetting consensus matrix. Was a bad letter provided?")
            }
            ans <- mat[row, ] != 0
            return(ans)
          })
          
#' @rdname sitesWithUnknown
setMethod("sitesWithUnknown",
          signature("DNAStringSet"),
          function(object){
            return(sitesHave(object, 'N'))
          })

#' @rdname sitesWithKnown
setMethod("sitesWithKnown",
          signature("DNAStringSet"),
          function(object){
            return(!sitesHave(object, 'N'))
          })

#' @rdname subsetSites
setMethod("subsetSites",
          signature("DNAStringSet", "integer"),
          function(object, index){
            index <- rep.int(list(index), length(object))
            return(object[index])
          }
)

#' @rdname subsetSequences
setMethod("subsetSequences",
          signature("DNAStringSet", "integer"),
          function(object, index){
            return(object[index])
          }
)