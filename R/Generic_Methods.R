

setGeneric("sequenceLength",
           function(object){
             standardGeneric("sequenceLength")
           }
)

setGeneric("statesPerBase",
          function(object) {
            standardGeneric("statesPerBase")
          }
)

setGeneric("polymorphicSites",
          function(object) {
            standardGeneric("polymorphicSites")
          }
)

setGeneric("sitesWithUnknown",
          function(object) {
            standardGeneric("sitesWithUnknown")
          }
)

setGeneric("sitesNotUnknown",
           function(object){
             standardGeneric("sitesNotUnknown")
           }
)

setGeneric("subsetSites",
          function(object, index){
            standardGeneric("subsetSites")
          }
)

setGeneric("subsetSequences",
           function(object, index){
             standardGeneric("subsetSequences")
           }
)

setGeneric("excludeUnknownSites",
           function(object) {
             standardGeneric("excludeUnknownSites")
           }
)

