### R code from vignette source 'Trees.Rnw'

###################################################
### code chunk number 1: Trees.Rnw:47-49
###################################################
options(width=70)
foo <- packageDescription("phangorn")


###################################################
### code chunk number 2: Trees.Rnw:65-67
###################################################
library(phangorn)
primates = read.phyDat("primates.dna", format="phylip", type="DNA")


###################################################
### code chunk number 3: Trees.Rnw:74-77
###################################################
dm = dist.dna(as.DNAbin(primates))
treeUPGMA = upgma(dm)
treeNJ = NJ(dm)


###################################################
### code chunk number 4: plotNJ
###################################################
layout(matrix(c(1,2), 2, 1), height=c(1,2))
par(mar = c(.1,.1,.1,.1))
plot(treeUPGMA, main="UPGMA")
plot(treeNJ, "unrooted", main="NJ")


###################################################
### code chunk number 5: figNJ
###################################################
getOption("SweaveHooks")[["fig"]]()
layout(matrix(c(1,2), 2, 1), height=c(1,2))
par(mar = c(.1,.1,.1,.1))
plot(treeUPGMA, main="UPGMA")
plot(treeNJ, "unrooted", main="NJ")


###################################################
### code chunk number 6: Trees.Rnw:99-101
###################################################
parsimony(treeUPGMA, primates)
parsimony(treeNJ, primates)


###################################################
### code chunk number 7: Trees.Rnw:104-107
###################################################
treePars = optim.parsimony(treeUPGMA, primates)
treeRatchet = pratchet(primates, trace = 0)
parsimony(c(treePars, treeRatchet), primates)


###################################################
### code chunk number 8: Trees.Rnw:110-111 (eval = FALSE)
###################################################
## (trees <- bab(subset(primates,1:10)))


###################################################
### code chunk number 9: Trees.Rnw:117-119
###################################################
fit = pml(treeNJ, data=primates)
fit


###################################################
### code chunk number 10: Trees.Rnw:122-123
###################################################
methods(class="pml")


###################################################
### code chunk number 11: Trees.Rnw:126-128
###################################################
fitJC = optim.pml(fit, TRUE)
logLik(fitJC)


###################################################
### code chunk number 12: Trees.Rnw:131-135
###################################################
fitGTR = update(fit, k=4, inv=0.2) 
fitGTR = optim.pml(fitGTR, TRUE,TRUE, TRUE, TRUE, TRUE, 
    control = pml.control(trace = 0))
fitGTR 


###################################################
### code chunk number 13: Trees.Rnw:138-139
###################################################
anova(fitJC, fitGTR) 


###################################################
### code chunk number 14: Trees.Rnw:142-144
###################################################
AIC(fitGTR) 
AIC(fitJC)


###################################################
### code chunk number 15: Trees.Rnw:147-148
###################################################
SH.test(fitGTR, fitJC) 


###################################################
### code chunk number 16: Trees.Rnw:151-152
###################################################
load("Trees.RData")


###################################################
### code chunk number 17: Trees.Rnw:154-155 (eval = FALSE)
###################################################
## mt = modelTest(primates)


###################################################
### code chunk number 18: Trees.Rnw:159-161
###################################################
library(xtable)
xtable(mt, caption="Summary table of modelTest", label="tab:modelTest")


###################################################
### code chunk number 19: Trees.Rnw:165-168
###################################################
env <- attr(mt, "env")
ls(envir=env)
(fit <- eval(get("HKY+G+I", env), env))


###################################################
### code chunk number 20: Trees.Rnw:172-174 (eval = FALSE)
###################################################
## bs = bootstrap.pml(fitJC, bs=100, optNni=TRUE, 
##     control = pml.control(trace = 0))


###################################################
### code chunk number 21: plotBS
###################################################
par(mar=c(.1,.1,.1,.1))
plotBS(fitJC$tree, bs)


###################################################
### code chunk number 22: figBS
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mar=c(.1,.1,.1,.1))
plotBS(fitJC$tree, bs)


###################################################
### code chunk number 23: Trees.Rnw:198-200
###################################################
options(prompt=" ")
options(continue="  ")


###################################################
### code chunk number 24: Trees.Rnw:202-225 (eval = FALSE)
###################################################
## library(parallel) # supports parallel computing
## library(phangorn)
## file="myfile"
## dat = read.phyDat(file)
## dm = dist.ml(dat)
## tree = NJ(dm)
## # as alternative for a starting tree:
## tree <- pratchet(dat) 
## 
## # 1. alternative: estimate an GTR model
## fitStart = pml(tree, dat, k=4, inv=.2)
## fit = optim.pml(fitStart, TRUE, TRUE, TRUE, TRUE, TRUE) 
##  
## # 2. alternative: modelTest  
## (mt <- modelTest(dat, multicore=TRUE)) 
## mt$Model[which.min(mt$BIC)]
## # choose best model from the table, assume now GTR+G+I
## env = attr(mt, "env")
## fitStart = eval(get("GTR+G+I", env), env) 
## fitStart = eval(get(mt$Model[which.min(mt$BIC)], env), env) 
## fit = optim.pml(fitStart, optNni=TRUE, optGamma=TRUE, optInv=TRUE, 
##     model="GTR")
## bs = bootstrap.pml(fit, bs=100, optNni=TRUE, multicore=TRUE)


###################################################
### code chunk number 25: Trees.Rnw:229-243 (eval = FALSE)
###################################################
## library(parallel) # supports parallel computing
## library(phangorn)
## file="myfile"
## dat = read.phyDat(file, type = "AA")
## dm = dist.ml(dat, model="JTT")
## tree = NJ(dm)
## 
## (mt <- modelTest(dat, model=c("JTT", "LG", "WAG"), multicore=TRUE)) 
## fitStart = eval(get(mt$Model[which.min(mt$BIC)], env), env) 
## 
## fitNJ = pml(tree, dat, model="JTT", k=4, inv=.2)
## fit = optim.pml(fitNJ, optNni=TRUE, optInv=TRUE, optGamma=TRUE)
## fit
## bs = bootstrap.pml(fit, bs=100, optNni=TRUE, multicore=TRUE)


###################################################
### code chunk number 26: Trees.Rnw:251-252
###################################################
toLatex(sessionInfo())


