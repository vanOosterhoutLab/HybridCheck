source("Simulation_Functions.R")

allResults <- lapply(seq(from = 1000, to = 10000, by = 1000),
                     function(x){
                       lapply(seq(from = 0, to = 10000, by = 1000),
                              function(y){
                                lapply(1:100, function(i) testATriplet(10e-6, x, y, TRUE))
                              })
                     })
save(allResults, file="SimulatedWithKnown.RData")
allResults2 <- lapply(allResults, function(x) lapply(x, function(y) do.call(rbind, y)))
allResults3 <- lapply(allResults2, function(x) do.call(rbind,x))
for(i in 1:length(allResults3)){
  allResults3[[i]]$ActualAgeBefore <- i * 1000
}
allResults4 <- do.call(rbind, allResults3)

# Convert Age (in number of generations) to mu * t for publication. 
allResults4$fiveAge <- allResults4$fiveAge * 10e-6
allResults4$fiftyAge <- allResults4$fiftyAge * 10e-6
allResults4$ninetyFiveAge <- allResults4$ninetyFiveAge * 10e-6
allResults4$ActualAge <- allResults4$ActualAge * 10e-6
allResults4$ActualAgeBefore <- allResults4$ActualAgeBefore * 10e-6

write.csv(allResults4, file = "KnownDatedAll.csv")

fiveAgeSummary <- ddply(allResults4, .(ActualAgeBefore, ActualAge), summarize,
                        mean = mean(fiveAge), sd = sd(fiveAge))
colnames(fiveAgeSummary) <- c("ActualAgeBefore", "ActualAge", "mean05", "sd05")
fiftyAgeSummary <- ddply(allResults4, .(ActualAgeBefore, ActualAge), summarize,
                         mean = mean(fiftyAge), sd = sd(fiftyAge))
colnames(fiftyAgeSummary) <- c("ActualAgeBefore", "ActualAge", "mean50", "sd50")
ninetyFiveAgeSummary <- ddply(allResults4, .(ActualAgeBefore, ActualAge), summarize,
                              mean = mean(ninetyFiveAge), sd = sd(ninetyFiveAge))
colnames(ninetyFiveAgeSummary) <- c("ActualAgeBefore", "ActualAge", "mean95", "sd95")
completeSummary <- cbind(fiveAgeSummary, fiftyAgeSummary[,3:4], ninetyFiveAgeSummary[,3:4])
write.csv(completeSummary, file = "KnownDatedSummary.csv")

ind <- which((completeSummary$ActualAgeBefore == 0.02) | (completeSummary$ActualAgeBefore == 0.04) | (completeSummary$ActualAgeBefore == 0.06) | (completeSummary$ActualAgeBefore == 0.08) | (completeSummary$ActualAgeBefore == 0.1))
summaryPlot <- ggplot(completeSummary[ind,], aes(x = ActualAge, y = mean50, shape = as.factor(ActualAgeBefore))) + 
  geom_point() + 
  geom_line() +
  geom_pointrange(aes(ymax = mean05, ymin = mean95))
