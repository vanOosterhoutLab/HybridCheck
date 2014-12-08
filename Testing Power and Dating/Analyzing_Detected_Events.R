source("Simulation_Functions.R")

allResultsDetection <- lapply(seq(from = 1000, to = 10000, by = 1000),
                              function(x){
                                lapply(seq(from = 0, to = 10000, by = 1000),
                                       function(y){
                                         lapply(1:100, function(i) testATriplet(10e-6, x, y, FALSE))
                                       })
                              })
save(allResultsDetection, file="SimulatedWithDetection.RData")
processedResultsDetection <- datesProcess(allResultsDetection)

# Convert Age (in number of generations) to mu * t for publication. 
processedResultsDetection$fiveAge <- processedResultsDetection$fiveAge * 10e-6
processedResultsDetection$fiftyAge <- processedResultsDetection$fiftyAge * 10e-6
processedResultsDetection$ninetyFiveAge <- processedResultsDetection$ninetyFiveAge * 10e-6
processedResultsDetection$ActualAge <- processedResultsDetection$ActualAge * 10e-6
processedResultsDetection$ActualAgeBefore <- processedResultsDetection$ActualAgeBefore * 10e-6

write.csv(processedResultsDetection, file = "DetectionDatedAll.csv")
detectionCompleteSummary <- summarizeDates(processedResultsDetection)
write.csv(detectionCompleteSummary, file = "DetectionDatedSummary.csv")
summaryPlots <- plotSummaryDates(detectionCompleteSummary)
ggsave("detectionDatesTotalSummary.png", plot = summaryPlots[[1]], width = 7,
       height = 7)

detectionsAsPercentage <- function(x, before, after){
  if(nrow(x) != 0){
    seqs <- unlist(strsplit(unique(as.character(x$Triplet)), ":"))
    correctpicks <- as.character(x$Sequence_Pair) == paste(list(seqs[1], seqs[2]), collapse = ":") | as.character(x$Sequence_Pair) == paste(list(seqs[1], seqs[3]), collapse = ":")
    x2 <- x[correctpicks,]
    if(nrow(x2) != 0){
      result <- data.frame(percentageDetected = sum((x2$Approximate_Length_BP)/50000)*100, AgeBefore = before, AgeAfter = after)
    } else {
      result <- data.frame(percentageDetected = 0, AgeBefore = before, AgeAfter = after)
    }
  } else {
    result <- data.frame(percentageDetected = 0, AgeBefore = before, AgeAfter = after)
  }
  return(result)
}

bootstrapCI <- function(pool, reps, upper, lower){
  bstrap <- c()
  for(b in 1:reps){
    bsamples <- sample(pool, length(pool), replace=TRUE)
    mbsamples <- mean(as.numeric(bsamples))
    bstrap <- c(bstrap, mbsamples)
  }
  lowerb <- quantile(bstrap, lower)
  upperb <- quantile(bstrap, upper)
  return(c(Lower = lowerb, Upper = upperb))
}

detectionPercentages <- lapply(1:length(allResultsDetection), function(x){
  lapply(1:length(allResultsDetection[[x]]), function(y){
    lapply(allResultsDetection[[x]][[y]], function(z){
      detectionsAsPercentage(z, x * 1000, (y-1) * 1000)
    })
  })
})
detectionPercentages2 <- lapply(detectionPercentages, function(x){
  lapply(x, function(z) do.call(rbind, z))
})
detectionPercentages3 <- lapply(detectionPercentages2, function(x) do.call(rbind, x))
detectionPercentagesTable <- do.call(rbind, detectionPercentages3)
detectionPercentagesTable$DivBefore <- detectionPercentagesTable$AgeBefore * 10e-6
detectionPercentagesTable$DivAfter <- detectionPercentagesTable$AgeAfter * 10e-6
write.csv(detectionPercentagesTable[,-c(2,3)], file = "detectionPercentagesFull.csv")
detectionSummaryTable <- ddply(detectionPercentagesTable, .(DivBefore, DivAfter), summarize,
                               mean = mean(percentageDetected), lower = bootstrapCI(percentageDetected, 100, 0.95, 0.05)[1], upper = bootstrapCI(percentageDetected, 100, 0.95, 0.05)[2])
write.csv(detectionSummaryTable, file = "detectionPercentagesSummary.csv")
ind <- which((detectionSummaryTable$DivBefore == 0.02) | (detectionSummaryTable$DivBefore == 0.04) | (detectionSummaryTable$DivBefore == 0.06) | (detectionSummaryTable$DivBefore == 0.08) | (detectionSummaryTable$DivBefore == 0.10))
detectionPlot <- ggplot(detectionSummaryTable[ind,], aes(x = DivAfter, y = mean, shape = as.factor(DivBefore))) +
  geom_point(size = 3) + geom_line() + geom_pointrange(aes(ymax = upper, ymin = lower)) +
  xlab(expression(paste("Age of the recombination event as ", mu*t))) +
  ylab("Percentage of recombinant regions correctly detected") +
  guides(shape = guide_legend(title = expression(paste("Divergence of parental sequences (", mu*t, ")"))))

ggsave("detectionPercentagePlot.png", plot = detectionPlot, width = 7, height = 5)
