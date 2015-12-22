# Classes for signal analyses:

# Generics and base classes.

#' Basic S4 class for storing the output and settings of a Recombination Signal Analysis.
#' 
#' @slot data A multiple sequence alignment stored as a DNAMultipleAlignment.
#' @slot results A RangedData variable that stores the recombination signal calculated
#' between sequences across the alignment, the location of said signal represented by IRanges.
setClass(
  "SignalAnalysis",
  representation(data = "DNAMultipleAlignment",
                 results = "RangedData")
)

#' A Virtual class to which all methods of computing a signal analysis should belong.
#' 
#' As a virtual class, has no fields, and should not be instantiated.
setClass("SignalAnalysisMethod",
         representation(),
         contains = "VIRTUAL")


#' Generic function used to perform a SignalAnalysis.
#' 
#' 
setGeneric("analyzeSignal", function(data, method) {
  standardGeneric("analyzeSignal")
})


# More specific methods of analysis.
setClass(
  "SlidingWindow",
  representation(width = "integer",
                 step = "integer"),
  contains = "SignalAnalysisMethod",
  prototype = prototype(width = 1000L,
                        step = 1000L)
)

setClass(
  "ReferenceWindowScan",
  representation(reference = "character"),
  contains = "SlidingWindow",
  prototype = prototype(width = 1000L,
                        step = 1000L)
)

setClass(
  "SSToReferenceWindowScan",
  representation = representation(settings = "ReferenceWindowScan"),
  contains = "SignalAnalysis"
)

setMethod("analyzeSignal",
          signature(data = "DNAMultipleAlignment", method = "ReferenceWindowScan"),
          function(data, method) {
            polymorphicOnly <- maskConservedSites(data, "union")
            results <-
              foreach(
                x = pairsRef(polymorphicOnly, method@reference), .combine = c, .multicombine = TRUE
              ) %dopar% {
                table <- slidingPoly(x, method@width, method@step)
                end(table) <-
                  unmaskedSites(x)[end(table)]
                start(table) <-
                  unmaskedSites(x)[start(table)]
                table$firstSeq <-
                  as.factor(method@reference)
                snames <- getSeqNames(x)
                table$secondSeq <-
                  as.factor(snames[snames != method@reference])
                return(table)
              }
            return(new(
              "SSToReferenceWindowScan",
              data = data,
              settings = method,
              results = results
            ))
          })

setMethod("sequenceLength",
          signature(object = "SignalAnalysis"),
          function(object){
            sequenceLength(object@data)
          })


