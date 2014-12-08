## ----eval=FALSE----------------------------------------------------------
#  vignette("HybRIDS_user_manual")

## ----eval=FALSE----------------------------------------------------------
#  install.packages("devtools", dep=T)
#  library(devtools)
#  install_github("Ward9250/HybRIDS", build_vignettes=TRUE)

## ----eval=FALSE----------------------------------------------------------
#  install_github("Ward9250/HybRIDS", ref="devel")

## ----message=FALSE-------------------------------------------------------
library(HybRIDS)

## ----eval=FALSE----------------------------------------------------------
#  MyAnalysis <- HybRIDS$new("~/Dropbox/MySequences.fas")

## ------------------------------------------------------------------------
data(MySequences)
MyAnalysis <- HybRIDS$new(MySequences)

## ------------------------------------------------------------------------
MyAnalysis

## ------------------------------------------------------------------------
MyAnalysis$showParameters("TripletGeneration")

## ------------------------------------------------------------------------
MyAnalysis$setParameters("TripletGeneration", DistanceThreshold = 0.1)

## ------------------------------------------------------------------------
library(HybRIDS)
MyAnalysis <- HybRIDS$new("~/Dropbox/MySequences.fas")
MyAnalysis$showParameters("TripletGeneration")

## ------------------------------------------------------------------------
# Let's make our groups:

FirstGroup <- c("Seq1", "Seq2", "Seq3")
SecondGroup <- c("Seq4", "Seq5", "Seq6")
ThirdGroup <- c("Seq7", "Seq8", "Seq9", "Seq10")
myGroups <- list(FirstGroup, SecondGroup, ThirdGroup) # Put the groups in a list.

# Let's edit the Triplet Generation settings
MyAnalysis$setParameters("TripletGeneration", Groups = myGroups, PartitionStrictness = 1L)

# All done, now let's print the settings again:
MyAnalysis$showParameters("TripletGeneration")

## ------------------------------------------------------------------------
MyAnalysis$setParameters("TripletGeneration", Method = c(1L, 2L),
                         Groups = myGroups, PartitionStrictness = 1L)

## ------------------------------------------------------------------------
tripletsToScan <- list(c("Seq1", "Seq4", "Seq7"), c("Seq1", "Seq4", "Seq8")) 

MyAnalysis$analyzeSS(tripletsToScan)

## ----error=FALSE---------------------------------------------------------
# This triplet is not allowed, because earlier the TripletGeneration
# settings were set so as only triplets which try to find recombination regions between 
tripletsToScan <- list(c("Seq1", "Seq2", "Seq3"))
MyAnalysis$analyzeSS(tripletsToScan)

## ------------------------------------------------------------------------
tripletsToSearch <- list(c("Seq1", "Seq4", "Seq7"), c("Seq1", "Seq4", "Seq8"))
MyAnalysis$findBlocks(tripletsToSearch)

## ----error=FALSE---------------------------------------------------------
tripletsToSearch <- list(c("Seq1", "Seq2", "Seq3"))
MyAnalysis$findBlocks(tripletsToSearch)

## ----error=FALSE---------------------------------------------------------
tripletsToSearch <- list(c("Seq1", "Seq4", "Seq9"))
MyAnalysis$findBlocks(tripletsToSearch)

# A triplet must be scanned for recombination signal first
MyAnalysis$analyzeSS(tripletsToSearch)
MyAnalysis$findBlocks(tripletsToSearch)

## ------------------------------------------------------------------------
tripletsToTabulate <- list(c("Seq1", "Seq4", "Seq7"), c("Seq1", "Seq4", "Seq8"))
MyAnalysis$tabulateDetectedBlocks(tripletsToTabulate)

## ------------------------------------------------------------------------
tripletsToDate <- list(c("Seq1", "Seq4", "Seq7"), c("Seq1", "Seq4", "Seq8"))
MyAnalysis$dateBlocks(tripletsToDate)

## ------------------------------------------------------------------------
MyAnalysis$tabulateDetectedBlocks(tripletsToDate)

## ------------------------------------------------------------------------
MyAnalysis$plotTriplets(c("Seq1", "Seq4", "Seq7"))[[1]]

## ------------------------------------------------------------------------
MyAnalysis$addUserBlock(c("Seq1", "Seq4"), firstbp = 287463, lastbp = 295890)
MyAnalysis$tabulateUserBlocks()
MyAnalysis$dateUserBlocks()
MyAnalysis$tabulateUserBlocks()

## ------------------------------------------------------------------------
MyAnalysis$clearUserBlocks(c("Seq1", "Seq4"))
MyAnalysis$tabulateUserBlocks()

