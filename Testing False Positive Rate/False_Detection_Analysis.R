# False Detection Analysis.
library(HybRIDS)

# First defining functions for analysis.

makeTriplets <- function(num=100, fileName, subpopSize){
  sequences <- read.dna(fileName, format="fasta")
  firsts <- sample(seq(from = 1, to = subpopSize), 100, replace = T)
  seconds <- sample(seq(from = subpopSize+1, to = subpopSize*2), 100, replace = T)
  thirds <- sample(seq(from = (subpopSize*2)+1, to = subpopSize*3), 100, replace = T)
  triplets <- lapply(1:length(firsts), function(i) sequences[c(firsts[i], seconds[i], thirds[i]),])
  seqnames <- lapply(triplets, function(x) labels(x))
  namechecks <- unlist(lapply(seqnames, function(x){grepl("_P1_", x[1]) && grepl("_P2_", x[2]) && grepl("_P3_", x[3])})) 
  if(!all(namechecks)){
    warning(paste("WARNING: Sequences have not been sampled to form triplets correctly from simuPOP output!!!\nTriplet: ", which(namechecks == FALSE)))
  }
  return(triplets)
}

HybRIDSanalysis <- function(triplet){
  HybRIDSobject <- HybRIDS$new(triplet)
  if(nrow(HybRIDSobject$DNA$InformativeSequence) < 3){
    return("Not enough unique sequences to form triplets")
  } else {
    HybRIDSobject$setParameters(Step="BlockDating", PValue = 0.05, BonfCorrection = TRUE,
                                MutationRate = 10e-6, MutationCorrection="JC69")
    HybRIDSobject$setParameters(Step="BlockDetection", AutoThresholds = TRUE, SDstringency = 2)
    HybRIDSobject$analyzeSS()
    HybRIDSobject$findBlocks()
    HybRIDSobject$dateBlocks()
    return(HybRIDSobject$tabulateDetectedBlocks())
  }
}

detectionsAsPercentage <- function(x){
  (sum(x$Approximate_Length_BP)/50000)*100
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

message("Running analysis to assess false detections, procedure scripted in False_Detection_Analysis.R")

# Make the sampled sets of three sequences from each DNA sequence file...
Founding1000Triplets <- makeTriplets(fileName = "FASTA Sequences/Founding1000PloidOne.fas", subpopSize = 500)
Founding2000Triplets <- makeTriplets(fileName = "FASTA Sequences/Founding2000PloidOne.fas", subpopSize = 500)
Founding3000Triplets <- makeTriplets(fileName = "FASTA Sequences/Founding3000PloidOne.fas", subpopSize = 500)
Founding4000Triplets <- makeTriplets(fileName = "FASTA Sequences/Founding4000PloidOne.fas", subpopSize = 500)
Founding5000Triplets <- makeTriplets(fileName = "FASTA Sequences/Founding5000PloidOne.fas", subpopSize = 500)
Founding6000Triplets <- makeTriplets(fileName = "FASTA Sequences/Founding6000PloidOne.fas", subpopSize = 500)
Founding7000Triplets <- makeTriplets(fileName = "FASTA Sequences/Founding7000PloidOne.fas", subpopSize = 500)
Founding8000Triplets <- makeTriplets(fileName = "FASTA Sequences/Founding8000PloidOne.fas", subpopSize = 500)
Founding9000Triplets <- makeTriplets(fileName = "FASTA Sequences/Founding9000PloidOne.fas", subpopSize = 500)
Founding10000Triplets <- makeTriplets(fileName = "FASTA Sequences/Founding10000PloidOne.fas", subpopSize = 500)

# Do HybRIDS analysis for many triplets...
resultsFounding1000 <- lapply(Founding1000Triplets, function(x) {HybRIDSanalysis(x)})
resultsFounding2000 <- lapply(Founding2000Triplets, function(x) {HybRIDSanalysis(x)})
resultsFounding3000 <- lapply(Founding3000Triplets, function(x) {HybRIDSanalysis(x)})
resultsFounding4000 <- lapply(Founding4000Triplets, function(x) {HybRIDSanalysis(x)})
resultsFounding5000 <- lapply(Founding5000Triplets, function(x) {HybRIDSanalysis(x)})
resultsFounding6000 <- lapply(Founding6000Triplets, function(x) {HybRIDSanalysis(x)})
resultsFounding7000 <- lapply(Founding7000Triplets, function(x) {HybRIDSanalysis(x)})
resultsFounding8000 <- lapply(Founding8000Triplets, function(x) {HybRIDSanalysis(x)})
resultsFounding9000 <- lapply(Founding9000Triplets, function(x) {HybRIDSanalysis(x)})
resultsFounding10000 <- lapply(Founding10000Triplets, function(x) {HybRIDSanalysis(x)})

falseDetectionsFounding1000 <- cbind(unlist(lapply(resultsFounding1000,
                                                   detectionsAsPercentage)), 1:100, 1000)
falseDetectionsFounding2000 <- cbind(unlist(lapply(resultsFounding2000,
                                                   detectionsAsPercentage)), 1:100, 2000)
falseDetectionsFounding3000 <- cbind(unlist(lapply(resultsFounding3000,
                                                   detectionsAsPercentage)), 1:100, 3000)
falseDetectionsFounding4000 <- cbind(unlist(lapply(resultsFounding4000,
                                                   detectionsAsPercentage)), 1:100, 4000)
falseDetectionsFounding5000 <- cbind(unlist(lapply(resultsFounding5000,
                                                   detectionsAsPercentage)), 1:100, 5000)
falseDetectionsFounding6000 <- cbind(unlist(lapply(resultsFounding6000,
                                                   detectionsAsPercentage)), 1:100, 6000)
falseDetectionsFounding7000 <- cbind(unlist(lapply(resultsFounding7000,
                                                   detectionsAsPercentage)), 1:100, 7000)
falseDetectionsFounding8000 <- cbind(unlist(lapply(resultsFounding8000,
                                                   detectionsAsPercentage)), 1:100, 8000)
falseDetectionsFounding9000 <- cbind(unlist(lapply(resultsFounding9000,
                                                   detectionsAsPercentage)), 1:100, 9000)
falseDetectionsFounding10000 <- cbind(unlist(lapply(resultsFounding10000,
                                                    detectionsAsPercentage)), 1:100, 10000)
falseDetectionsResults <- rbind(falseDetectionsFounding1000, falseDetectionsFounding2000,
                                falseDetectionsFounding3000, falseDetectionsFounding4000,
                                falseDetectionsFounding5000, falseDetectionsFounding6000,
                                falseDetectionsFounding7000, falseDetectionsFounding8000,
                                falseDetectionsFounding9000, falseDetectionsFounding10000)
colnames(falseDetectionsResults) <- c("PercentageFalselyDetected", "Triplet", "Divergence")
falseDetectionsResults <- as.data.frame(falseDetectionsResults)
falseDetectionsResults$Divergence <- falseDetectionsResults$Divergence * 10e-6

falseDetectionPlot <- ggplot(falseDetectionsResults, aes(x=Divergence, y=PercentageFalselyDetected)) + geom_point()
ggsave(paste(ResultsDirectory, "FalseDetectionsPlot.png", sep=""), falseDetectionPlot)
ggsave(paste(ResultsDirectory, "FalseDetectionsPlot.pdf", sep=""), falseDetectionPlot)

# Summarize the data points...
sets <- c("falseDetectionsFounding1000", "falseDetectionsFounding2000",
          "falseDetectionsFounding3000", "falseDetectionsFounding4000",
          "falseDetectionsFounding5000", "falseDetectionsFounding6000",
          "falseDetectionsFounding7000", "falseDetectionsFounding8000",
          "falseDetectionsFounding9000", "falseDetectionsFounding10000")

falseDetectionMeans <- cbind(
  unlist(lapply(sets, function(x){mean(get(x)[,1])})),
  c(1000,2000,3000,4000,5000,6000,7000,8000,9000,10000) * 10e-6)

falseDetectionCIs <- do.call(rbind,
                             lapply(sets, function(x){
                               bootstrapCI(get(x)[,1], 1000, .975, .025)
                               }))

Xlab1 <- expression(paste("Divergence (", mu*t, ") between populations."))
XaxisScale <- scale_x_continuous(breaks = seq(from=1000, to=50000, by=1000), labels = seq(from=1000*10e-6, to=50000*10e-6, by=1000*10e-6))
YaxisScale <- scale_y_continuous(limits=c(0,100), breaks = seq(from=0, to=100, by=10), labels = seq(from=0, to=1, by=0.1))
YaxisScale2 <- scale_y_continuous(breaks = seq(from=0, to=100, by=1), labels = seq(from=0, to=1, by=0.01))

falseDetectionSummary <- data.frame(cbind(falseDetectionMeans, falseDetectionCIs))
colnames(falseDetectionSummary) <- c("Mean", "Divergence", "lowerCI", "upperCI")
limits <- aes(ymax = upperCI, ymin = lowerCI)
falseSummaryPlot <- ggplot(falseDetectionSummary, aes(x=Divergence, y=Mean)) + geom_point() + geom_line() + 
  geom_errorbar(limits) + 
  ylab("Proportion of sequence incorrectly identified as recombination.") + 
  xlab(Xlab1) +
  XaxisScale
falseSummaryPlot2 <- falseSummaryPlot + YaxisScale
falseSummaryPlot3 <- falseSummaryPlot + YaxisScale2
ggsave(paste(ResultsDirectory, "Plots/FalseDetectionsSummary.png", sep=""), falseSummaryPlot3)
ggsave(paste(ResultsDirectory, "Plots/FalseDetectionsSummary.pdf", sep=""), falseSummaryPlot3)

write.csv(falseDetectionSummary, paste(ResultsDirectory, "summaryFalsePos.csv", sep=""))
write.csv(falseDetectionsResults, paste(ResultsDirectory, "fullFalsePos.csv", sep=""))