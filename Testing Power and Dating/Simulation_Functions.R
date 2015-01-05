#!/usr/bin/Rscript
library(HybRIDS)
library(phangorn)
library(plyr)

analyzeWithKnown <- function(triplet){
  hyb <- HybRIDS$new(triplet)
  if(hyb$DNA$numberOfSequences() < 3){
    return("Not enough unique sequences to form triplets")
  } else {
    hyb$setParameters(Step="BlockDating", PValue = 0.05, BonfCorrection = TRUE,
                      MutationRate = 10e-6, MutationCorrection="JC69", DateAnyway = TRUE)
    hyb$addUserBlock(paste0(hyb$DNA$getSequenceNames()[1], ":",
                            hyb$DNA$getSequenceNames()[2]), 1, 30000)
    hyb$addUserBlock(paste0(hyb$DNA$getSequenceNames()[1],
                            ":", hyb$DNA$getSequenceNames()[3]), 30001, 50000)
    hyb$dateUserBlocks()
    return(hyb$tabulateUserBlocks())
  }
}

analyzeWithDetection <- function(triplet){
  hyb <- HybRIDS$new(triplet)
  if(hyb$DNA$numberOfSequences() < 3){
    return("Not enough unique sequences to form triplets")
  } else {
    hyb$setParameters(Step="BlockDating", PValue = 0.05, BonfCorrection = TRUE,
                      MutationRate = 10e-6, MutationCorrection="JC69")
    hyb$setParameters(Step="BlockDetection", AutoThresholds = TRUE, SDstringency = 2)
    hyb$analyzeSS()
    hyb$findBlocks()
    hyb$dateBlocks()
    return(hyb$tabulateDetectedBlocks())
  }
}

testATriplet <- function(mu, ngenbefore, ngenafter, knownBlocks){
  divBefore <- 2*mu*ngenbefore
  divAfter <- 2*mu*ngenafter
  divBeforeGuideTree <- compute.brlen(stree(2), divBefore/2)
  parentalSequences <- as.character(simSeq(divBeforeGuideTree, l = 50000))
  recombinantOffspring <- c(parentalSequences[1,1:30000],
                                      parentalSequences[2,30001:50000])
  triplet <- rbind(parentalSequences, recombinantOffspring)
  divAfterGuideTree <- compute.brlen(stree(2), divAfter/2)
  agedRecombinant <- as.character(simSeq(divAfterGuideTree, l = 50000, rootseq = triplet[3,]))[1,]
  agedParent1 <- as.character(simSeq(divAfterGuideTree, l = 50000, rootseq = triplet[1,]))[1,]
  agedParent2 <- as.character(simSeq(divAfterGuideTree, l = 50000, rootseq = triplet[2,]))[1,]
  tripletToAnalyze <- as.DNAbin(rbind(agedRecombinant, agedParent1, agedParent2))
  statsTable <- data.frame(Mu = mu,
                           GenBefore = ngenbefore,
                           setDivBefore = divBefore,
                           divBefore = dist.dna(as.DNAbin(parentalSequences), model="JC69")[1],
                           GenAfter = ngenafter,
                           setDivAfter = divAfter,
                           parentalDivAfter = dist.dna(tripletToAnalyze, model="JC69")[3],
                           recombinantParent1Div = dist.dna(tripletToAnalyze, model="JC69")[1],
                           recombinantParent2Div = dist.dna(tripletToAnalyze, model="JC69")[2],
                           rp1DivInRegion = dist.dna(tripletToAnalyze[,1:30000], model="JC69")[1],
                           rp2DivInRegion = dist.dna(tripletToAnalyze[,30001:50000], model="JC69")[2])
  message(paste0('Simulation Summary\n==================\n\n', 'Mutation Rate: ', mu,
                 '\n\nEvolution between parents Before Recombination:\n\tNumber of Generations: ',
                 statsTable$GenBefore,
                 '\n\tAmount of Divergence between parents set: ', statsTable$setDivBefore,
                 '\n\tActual Divergence between parents: ',
                 statsTable$divBefore,
                 '\n\nAmount of Evolution after recombination:\nNumber of Generations: ',
                 statsTable$GenAfter, '\nAmount of Divergence between sequences set: ', statsTable$setDivAfter,
                 '\nActual Divergence between parents: ',
                 statsTable$parentalDivAfter,
                 '\nActual Divergence between recombinant and aged parent 1: ',
                 statsTable$recombinantParent1Div,
                 '\nActual Divergence between recombinant and aged parent 2: ',
                 statsTable$recombinantParent2Div,
                 '\nActual Divergence between recombinant and aged parent 1',
                 ' in the recombinant region between them: ',
                 statsTable$rp1DivInRegion,
                 '\nActual Divergence between recombinant and aged parent 2',
                 ' in the recombinant region between them: ',
                 statsTable$rp1DivInRegion))
  message("Now analyzing the simulated triplet with HybRIDS...")
  if(knownBlocks){
    hybridsResults <- analyzeWithKnown(tripletToAnalyze)
  } else {
    hybridsResults <- analyzeWithDetection(tripletToAnalyze)
  }
  if(nrow(hybridsResults) > 0){
    hybridsResults$ActualAge <- ngenafter
    hybridsResults$AgeBefore <- ngenbefore
  }
  return(hybridsResults)     
}

plotSummaryDates <- function(summaryDf){
  yLabel <- expression(paste("Estimated age of recombinant regions expressed as ", mu*t))
  xLabel <- expression(paste("Actual age of recombinant regions expressed as ", mu*t))
  Ylab <- ylab(yLabel)
  Xlab <- xlab(xLabel)
  Legend <- guides(colour=guide_legend(title="Parental sequence divergence before recombination"))
  summaryPlot <- ggplot(summaryDf, 
                        aes(x = ActualAge, y = mean50, colour = as.factor(ActualAgeBefore))) +
    geom_point() + geom_line() + geom_pointrange(aes(ymax = mean05, ymin = mean95)) +
    geom_abline(intercept = 0, slope = 1) + Ylab + Xlab +
    Legend
  summaryPlot5 <- ggplot(summaryDf, 
                         aes(x = ActualAge, y = mean05, colour = as.factor(ActualAgeBefore))) +
    geom_point() + geom_line() + geom_pointrange(aes(ymax = mean05 + sd05, ymin = mean05 - sd05)) +
    geom_abline(intercept = 0, slope = 1) + Ylab + Xlab +
    Legend
  summaryPlot50 <- ggplot(summaryDf, 
                          aes(x = ActualAge, y = mean50, colour = as.factor(ActualAgeBefore))) +
    geom_point() + geom_line() + geom_pointrange(aes(ymax = mean50 + sd50, ymin = mean50 - sd50)) +
    geom_abline(intercept = 0, slope = 1) + Ylab + Xlab +
    Legend
  summaryPlot95 <- ggplot(summaryDf, 
                          aes(x = ActualAge, y = mean95, colour = as.factor(ActualAgeBefore))) +
    geom_point() + geom_line() + geom_pointrange(aes(ymax = mean95 + sd95, ymin = mean95 - sd95)) +
    geom_abline(intercept = 0, slope = 1) + Ylab + Xlab + Legend
  return(list(totalSummary = summaryPlot, summary05 = summaryPlot5,
              summary50 = summaryPlot50, summary95 = summaryPlot95))
}

datesProcess <- function(simulationOutput){
  simulationOutput2 <- lapply(simulationOutput, function(x) lapply(x, function(y) do.call(rbind, y)))
  simulationOutput3 <- lapply(simulationOutput2, function(x) do.call(rbind,x))
  for(i in 1:length(simulationOutput3)){
    simulationOutput3[[i]]$ActualAgeBefore <- i * 1000
  }
  simulationOutput4 <- do.call(rbind, simulationOutput3)
  return(simulationOutput4)
}

summarizeDates <- function(dataFrame){
  fiveAgeSummary <- ddply(dataFrame, .(ActualAgeBefore, ActualAge), summarize,
                          mean = mean(fiveAge), sd = sd(fiveAge))
  colnames(fiveAgeSummary) <- c("ActualAgeBefore", "ActualAge", "mean05", "sd05")
  fiftyAgeSummary <- ddply(dataFrame, .(ActualAgeBefore, ActualAge), summarize,
                           mean = mean(fiftyAge), sd = sd(fiftyAge))
  colnames(fiftyAgeSummary) <- c("ActualAgeBefore", "ActualAge", "mean50", "sd50")
  ninetyFiveAgeSummary <- ddply(dataFrame, .(ActualAgeBefore, ActualAge), summarize,
                                mean = mean(ninetyFiveAge), sd = sd(ninetyFiveAge))
  colnames(ninetyFiveAgeSummary) <- c("ActualAgeBefore", "ActualAge", "mean95", "sd95")
  completeSummary <- cbind(fiveAgeSummary, fiftyAgeSummary[,3:4], ninetyFiveAgeSummary[,3:4])
  return(completeSummary)
}

