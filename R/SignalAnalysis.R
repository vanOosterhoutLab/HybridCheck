# Classes for signal analyses:

# Generics and base classes.

#' Basic S4 class for storing the output and settings of a Recombination
#' Signal Analysis.
#' 
#' @slot data A multiple sequence alignment stored as a DNAMultipleAlignment.
#' @slot results A RangedData variable that stores the recombination signal
#' calculated between sequences across the alignment, the location of said 
#' signal represented by IRanges.
setClass(
  "SignalAnalysis",
  representation(data = "DNAMultipleAlignment",
                 results = "RangedData")
)


setMethod("sequenceLength",
          signature(object = "SignalAnalysis"),
          function(object) {
            sequenceLength(object@data)
          })


#' Generic function used to perform a SignalAnalysis.
#' 
#' Analyse the recombination signal present in a variable containing
#' sequence alignment data.
#' 
#' The specific method used to analyse the data is dispatched based on the
#' type of the input data, and the class of the method veriable.
#' 
#' @field data A variable containing the sequence data you wish to analyse.
#' @field method A variable which inherits from SignalAnalysisMethod. This
#' variable needs to contain all the options for the given SignalAnalysis
#' method.
#' @export
setGeneric("analyzeSignal", function(data, method) {
  standardGeneric("analyzeSignal")
})

#' A Virtual class to which all methods of computing a signal analysis should
#' belong.
#' 
#' As a virtual class, has no fields, and should not be instantiated.
setClass("SignalAnalysisMethod",
         representation(),
         contains = "VIRTUAL")

#' A class for a Sliding Window based method of analysing recombination signal.
#' 
#' This class, and others which inherit from it, represent methods of analysing
#' recombination signal which utilise sliding windows to copute their values.
#' For example, the class ReferenceWindowScan inherits from this class.
#' @slot width An integer specifying the width of the sliding window in the 
#' analysis in base pairs.
#' @slot step An integer specifying the number of base positions by which the
#' sliding window in the analysis will move between each calculation.
setClass(
  "SlidingWindow",
  slots = c(width = "integer",
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
#' @slot width An integer specifying the width of the sliding window in the
#' analysis in base pairs.
#' @slot step An integer specifying the number of base positions by which the 
#' sliding window in the analysis will move between each calculation.
#' @export
setClass(
  "ReferenceWindowScan",
  slots = c(reference = "character"),
  contains = "SlidingWindow",
  prototype = prototype(width = 1000L,
                        step = 1000L)
)


#' Sliding window scan of recombination signal between many suery sequences 
#' and one reference sequence.
#' 
#' This class represents and analysis of recombination signal computed by 
#' calculating sequence similarity across sliding windows, between one
#' reference sequence, and many query sequences. As such this class
#' inherits from ReferenceWindowScan, which inherits from SlidingWindow.
#' @slot settings The ReferenceWindowScan object used to provide settings
#' for the performed analysis.
setClass(
  "SSToReferenceWindowScan",
  slots = c(settings = "ReferenceWindowScan"),
  contains = "SignalAnalysis"
)


#' @rdname analyzeSignal
setMethod("analyzeSignal",
          signature(data = "DNAMultipleAlignment", method = "ReferenceWindowScan"),
          function(data, method) {
            noNs <- maskNs(data, "union")
            polymorphicOnly <- maskConservedSites(noNs, "union")
            results <- foreach(x = pairsRef(polymorphicOnly, ref = method@reference, checkFunc = function(x) TRUE),
                               .combine = c,
                               .multicombine = TRUE) %dopar% {
              table <- slidingPoly(x, method@width, method@step)
              end(table) <- unmaskedSites(x)[end(table)]
              start(table) <- unmaskedSites(x)[start(table)]
              table$firstSeq <- as.factor(method@reference)
              snames <- getSeqNames(x)
              table$secondSeq <- as.factor(snames[snames != method@reference])
              return(table)
                               }
            return(new("SSToReferenceWindowScan", data = data, settings = method, results = results))
          })


#' @rdname analyzeSignal
setMethod("analyzeSignal",
          signature(data = "DNAStringSet", method = "ReferenceWindowScan"),
          function(data, method){
            known <- sitesWithKnown(data)
            polymorphic <- polymorphicSites(data)
            informativeIdx <- which(known & polymorphic)
            infSubset <- subsetSites(data, informativeIdx)
            itr <- pairsRef(infSubset, ref = method@reference)
            results <- foreach(x = itr,
                               .combine = c,
                               .multicombine = TRUE) %do% {
                                 table <- slidingPoly(x, method@width, method@step)
                                 end(table) <- informativeIdx[end(table)]
                                 start(table) <- informativeIdx[start(table)]
                                 return(table)
                               }
            
            return(l)
          }
)
