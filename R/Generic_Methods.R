
#' @name sequenceLength
#' @title Get the length of a sequence alignment
#' @description Generic method to get the sequence length from an object.
#' @rdname sequenceLength
#' @param object An object which represents sequence data.
#' @export
setGeneric("sequenceLength",
           function(object){
             standardGeneric("sequenceLength")
           }
)

#' @name statesPerBase
#' @title Compute the number of states at each position in sequence data.
#' @description Generic method to get the number of states at each position in objects containing sequence data.
#' @rdname statesPerBase
#' @param object An object which represents sequence data.
#' @export
setGeneric("statesPerBase",
          function(object) {
            standardGeneric("statesPerBase")
          }
)

#' @name polymorphicSites
#' @title Compute which sites in DNA sequence data are polymorphic.
#' @description Generic method to get the indicies of polymorphic sites in DNA sequence data objects.
#' @rdname polymorphicSites
#' @param object An object which represents sequence data.
#' @export
setGeneric("polymorphicSites",
          function(object) {
            standardGeneric("polymorphicSites")
          }
)

#' @name sitesHave
#' @title Compute which sites in DNA sequence data have a given letter.
#' @description Generic method to get the indicies of sites in DNA sequence data objects that have a certain letter.
#' @rdname sitesHave
#' @param object An object which represents sequence data.
#' @param letter The letter to look for at each site.
#' @export
setGeneric("sitesHave",
           function(object, letter){
             standardGeneric("sitesHave")
           }
)

#' @name sitesWithUnknown
#' @title Compute which sites in DNA sequence data have 'N's.
#' @description Generic method to get the indicies of unknown sites in DNA sequence data objects.
#' @rdname sitesWithUnknown
#' @param object An object which represents sequence data.
#' @export
setGeneric("sitesWithUnknown",
          function(object) {
            standardGeneric("sitesWithUnknown")
          }
)

#' @name sitesWithKnown
#' @title Compute which sites in DNA sequence data don't have 'N's.
#' @description Generic method to get the indicies of known sites in DNA sequence data objects.
#' @rdname sitesWithKnown
#' @param object An object which represents sequence data.
#' @export
setGeneric("sitesWithKnown",
           function(object){
             standardGeneric("sitesWithKnown")
           }
)

#' @name subsetSites
#' @title Generic method Create an DNA sequence object which only contains a subset of positions.
#' @description  Create an DNA sequence object which only contains the positions indicated by index.
#' @rdname subsetSites
#' @param object An object which represents sequence data.
#' @param index A numeric vector of sites to subset.
#' @export
setGeneric("subsetSites",
          function(object, index){
            standardGeneric("subsetSites")
          }
)

#' @name subsetSequences
#' @title Generic method Create an DNA sequence object which only contains a subset of sequences.
#' @description Create an DNA sequence object which only contains the sequences by index.
#' @rdname subsetSequences
#' @param object An object which represents sequence data.
#' @param index A numeric vector of sites to subset.
#' @export
setGeneric("subsetSequences",
           function(object, index){
             standardGeneric("subsetSequences")
           }
)

#' @name excludeUnknownSites
#' @title Subset sequence data objects such that they no longer contain N's
#' @description Generic method to eliminate sites with 'N's from sequence data.
#' @rdname excludeUnknownSites
#' @param object An object which represents sequence data.
#' @export
setGeneric("excludeUnknownSites",
           function(object) {
             standardGeneric("excludeUnknownSites")
           }
)

