# Signal Visualisation.

setClass("VisualisationMethod",
         contains = "VIRTUAL")

setClass("HeatPlot",
         representation = representation(tileScale = "integer"),
         prototype = prototype(tileScale = 100L),
         contains = "VisualisationMethod")

setClass("LinesPlot",
         contains = "VisualisationMethod")

setGeneric("visualize", function(obj, plotMethod) {
  standardGeneric("visualize")
})

setMethod("visualize",
          signature(obj = "SSToReferenceWindowScan", plotMethod = "HeatPlot"),
          function(obj, plotMethod) {
            plotSegments <-
              calculatePlotSlices(sequenceLength(obj), plotMethod@tileScale)
            plotvalues <-
              foreach(
                x = rangedDataSpaces(obj@results, space(obj@results)), .combine = cbind, .multicombine = TRUE
              ) %dopar% {
                return(calculateSliceValues(plotSegments, x))
              }
            colnames(plotvalues) <- unique(obj@results$secondSeq)
            plottingTable <- cbind(plotSegments, plotvalues)
            plottingTable <-
              melt(plottingTable, id.vars = colnames(plottingTable)[1:3])
            heatplot <-
              ggplot(plottingTable, aes(x = bpCenter, y = variable)) +
              geom_tile(aes(fill = value), colour = "black") +
              scale_fill_gradient(
                name = "Sequence Similarity %"
              ) +
              xlab("Base position") +
              ylab("Query Sequence")
            return(heatplot)
          })

setMethod("visualize",
          signature(obj = "SSToReferenceWindowScan", plotMethod = "LinesPlot"),
          function(obj, plotMethod) {
            plotTable <- data.frame(signal = obj@results$signal,
                                    bpCenter = mid(obj@results$ranges),
                                    querySeq = obj@results$secondSeq)
            lineplot <- ggplot(plotTable, aes(x = bpCenter, y = signal, colour = querySeq)) +
              geom_line() +
              xlab("Base Position") +
              ylab("Sequence Similarity %")
            return(lineplot)
          })

setGeneric("calculatePlotSlices", function(length, scale){
  standardGeneric("calculatePlotSlices")
})

setMethod("calculatePlotSlices",
          signature(length = "integer", scale = "integer"),
          function(length, scale){
            div <- length / scale
            frame <- data.frame(bpStart = seq(from = 1, to = length, by = div),
                                bpEnd = seq(from=div, to = length, by = div)) 
            frame$bpCenter <- round(frame$bpStart +  (div / 2))
            return(frame)
          }
)

setGeneric("calculateSliceValues", function(plotFrame, signal){
  standardGeneric("calculateSliceValues")
})

setMethod("calculateSliceValues",
          signature(plotFrame = "data.frame", signal = "RangedData"),
          function(plotFrame, signal) {
            rowit <- iapply(plotFrame, 1)
            values <- foreach(x = rowit, .combine = c) %do% {
              bool1 <- start(signal) <= x$bpEnd
              bool2 <- x$bpStart <= end(signal)
              index <- which(bool1 & bool2)
              return(mean(signal$signal[index]))
            }
            return(values)
          })








