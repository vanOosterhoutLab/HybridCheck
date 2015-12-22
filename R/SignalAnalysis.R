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


setMethod("sequenceLength",
          signature(object = "SignalAnalysis"),
          function(object){
            sequenceLength(object@data)
          })



#' A Virtual class to which all methods of computing a signal analysis should belong.
#' 
#' As a virtual class, has no fields, and should not be instantiated.
setClass("SignalAnalysisMethod",
         representation(),
         contains = "VIRTUAL")


#' Generic function used to perform a SignalAnalysis.
#' 
#' @field data A variable containing the sequence data you wish to analyse.
#' @field method A variable which inherits from SignalAnalysisMethod. This
#' variable needs to contain all the options for the given SignalAnalysis method.
#' The specific method analyzeSignal is then al 
setGeneric("analyzeSignal", function(data, method) {
  standardGeneric("analyzeSignal")
})


#' A class for a Sliding Window based method of analysing recombination signal.
#' 
#' This class, and others which inherit from it, represent methods of analysing
#' recombination signal which utilise sliding windows to copute their values.
#' For example, the class ReferenceWindowScan inherits from this class.
#' @slot width An integer specifying the width of the sliding window in the analysis
#' in base pairs.
#' @slot step An integer specifying the number of base positions by which the sliding window in the 
#' analysis will move between each calculation.
setClass(
  "SlidingWindow",
  representation(width = "integer",
                 step = "integer"),
  contains = "SignalAnalysisMethod",
  prototype = prototype(width = 1000L,
                        step = 1000L)
)


#' A class for a Sliding Window based method of analysing recombination signal,
#' by comparing sequence similarity of many sequences, against one reference.
#' 
#' This class, and others which inherit from it, represent methods of analysing
#' recombination signal which utilise sliding windows to copute their values.
#' For example, the class ReferenceWindowScan inherits from this class.
#' @slot width An integer specifying the width of the sliding window in the analysis
#' in base pairs.
#' @slot step An integer specifying the number of base positions by which the sliding window in the 
#' analysis will move between each calculation.
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


