#' A Reference Class for performing and storing results of ABBA-BABA tests.
#' @name ABBABABA
#' @field results A list of reference objects of type ABBABABArecord.
FTTester <- setRefClass("FTTester",
                         
                         fields = list(genomicOrRegion = "character",
                                       results = "list"),
                         
                         methods = list(
                           
                           doTest =
                             function(dna){
                               "Performs the ABBA-BABBA test, accepts a 4 sequence alignment as a Biostrings object."
                               if(length(dna) != 4){stop("ERROR: Alignment does not have 4 sequences.")}
                               observedAbbaBabaSites <- getAbbaBabaSites(dna)
                               
                               
                               
                               
                               
                               
                               numberOfAbba <- length(abbaSites)
                               numberOfBaba <- length(babaSites)
                               D <- (numberOfAbba - numberOfBaba) / (numberOfAbba + numberOfBaba)
                               
                               if(numberOfAbba < numberOfBaba){
                                 success <- numberOfAbba
                               } else {
                                 success <- numberOfBaba
                               }
                               pValue <- pbinom(success, (numberOfAbba + numberOfBaba), .5)
                               results <<- append(results, ABBABABArecord$new(      ))
                             },
                           
                           
                           
                           
                           )
                         )


getAbbaBabaSites <- function(dna){
  twoStates <- which(colSums(consensusMatrix(dna) != 0) == 2)
  twoStatesSeq <- DNAStringSet(character(length = 4))
  for(i in 1:4){
    twoStatesSeq[[i]] <- dna[[i]][twoStates]
  }
  oneNotTwo <- colSums(consensusMatrix(twoStatesSeq[c(1, 2)]) != 0) > 1
  threeNotFour <- colSums(consensusMatrix(twoStatesSeq[c(3, 4)]) != 0) > 1
  abbaBabbaSites <- which(oneNotTwo & threeNotFour)
  abbaBabbaSeq <- DNAStringSet(character(length = 4))
  for(i in 1:4){
    abbaBabbaSeq[[i]] <- twoStatesSeq[[i]][abbaBabbaSites]
  }
  abbaSites <- which(colSums(consensusMatrix(abbaBabbaSeq[c(2,3)]) != 0) == 1)
  babaSites <- which(colSums(consensusMatrix(abbaBabbaSeq[c(1,3)]) != 0) == 1)
  whichSitesAbba <- twoStates[abbaBabbaSites[abbaSites]]
  whichSitesBaba <- twoStates[abbaBabbaSites[babaSites]]
  return(list(AbbaSites = whichSitesAbba, BabaSites = whichSitesBaba))
}


FTTestRecord <- retRefClass("FTTestRecord",
                              
                              fields = list(
                                sequences = "character",
                                globalD = "numeric"
                                
                                
                                ),
                              
                              methods = list()
                              )




CalcD <- function(alignment = "alignment.fasta"){
  # Use seqinr's read.alignment function.
  alignment <- read.alignment(alignment, format = "fasta")
  
  # Convert the alignment to a matrix of characters.
  alignment.matrix <- matrix(, length(alignment$nam), nchar(alignment$seq[[1]]))
  for(i in 1:length(alignment$nam)){
    alignment.matrix[i, ] <- unlist(strsplit(alignment$seq[[i]], ""))
  }
  
  abba <- 0
  baba <- 0
  
  for(i in 1:ncol(alignment.matrix)){                                                    # For each site in the alignment....
    if(length(unique(alignment.matrix[, i])) == 2){                                        # IF there are precisely two different states...
      if(alignment.matrix[1, i] != alignment.matrix[2, i]){                                  # IF seq 1 DOES NOT EQUAL seq2 at this site... 
        if(alignment.matrix[4, i] != alignment.matrix[3, i]){                                  # IF seq 4 DOES NOT EQUAL seq3 at this site...
          if(alignment.matrix[3, i] == alignment.matrix[1, i]) {baba <- baba + 1}                # IF seq 3 DOES EQUAL seq1 at this site... IT IS BABA.
          if(alignment.matrix[2, i] == alignment.matrix[3, i]) {abba <- abba + 1}                # IF seq 2 DOES EQUAL seq3 at this site... IT IS ABBA.
        }
      }
    }
  }
  d <- (abba - baba) / (abba + baba) #what its all about
  if(abba < baba) success<-abba
  else success<-baba
  user.result <- list()
  user.result$pval <- pbinom(success,(abba+baba),.5,lower.tail=T)
  user.result$d.stat <- d
  user.result$align.length <- ncol(alignment.matrix) - 1
  user.result$abba.sites <- abba
  user.result$baba.sites <- baba
  print(paste("Sites in alignment =", ncol(alignment.matrix)))
  print(paste("Number of sites with ABBA pattern =", abba))
  print(paste("Number of sites with BABA pattern =", baba))
  print(paste("D statistic =", d))
  print(paste("P-value =", user.result$pval, " based on the binomial distribution"))
  return(user.result)
} 






abbaSeq <- DNAStringSet(character(length = 4))
for(i in 1:4){
  abbaSeq[[i]] <- dna[[i]][whichSitesAbba]
}

babaSeq <- DNAStringSet(character(length = 4))
for(i in 1:4){
  babaSeq[[i]] <- dna[[i]][whichSitesBaba]
}