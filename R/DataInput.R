# A function for reading DNA data into HybRIDS from FASTA format.
# Last edited by Ben. J. Ward on 02/04/2013. 

#' HybRIDS - A fast, simple to use, and attractive recombination detection and dating package.
#' @import ggplot2 grid gridExtra xtable png gWidgets gWidgetstcltk ape
#' @importFrom seqinr read.fasta
#' @docType package
#' @name HybRIDS
NULL


#'A function for loading in DNA sequences from fasta format and convert to HybRIDS format.
#'
#'The function for loading in and converting DNA sequence data in FASTA format. It can also perform a 
#'pre-analysis of whole sequence similarity to try and remove sequences unlikely to contain
#'recombinant regions e.g. those sequences that may be too similar. Note that analysis functions downstream of this function
#'are somewhat intelligent in that they know if there are no informative sites.
#'@export
dna.data.prepare <- function(x=NULL, method=1) {
  if(is.null(x) == TRUE){
    cat("
        Enter the correct filepath for your fasta format alignment file below.
        For example on a linux machine: '~/Bob/Documents/Data/sequences.fasta'.
        If the working directory of R is set to the folder in which your datafile is stored,
        then you only need to type the name, e.g. 'sequences.fasta'.
        ")
    dataFile <- readline("File Path: ")
  } else {
    dataFile <- x
  }
  dnaSeqLst <- read.fasta(file = dataFile, 
                          seqtype = "DNA", as.string = FALSE, forceDNAtolower = TRUE,
                          set.attributes = TRUE, legacy.mode = TRUE, seqonly = FALSE, strip.desc = FALSE)
  dnaSeqLengths <- sapply(dnaSeqLst, length)
  if(all(dnaSeqLengths == dnaSeqLengths[[1]]) == FALSE) stop("\nAll sequences are not of the same length. Please check your input files...\nInput files must be fasta files of an alignment, and 'not' the raw sequence files.\nAborting...")
  cat("\nAll sequences are of the same length - good - continuing process...")
  cat("\nFormatting data...")
  completeDNA <- do.call(rbind, dnaSeqLst)
  colnames(completeDNA) <- 1:dnaSeqLengths[[1]]
  rownames(completeDNA) <- names(dnaSeqLst)
  cat("\nDNA data formatted successfully!")
  combos <- combn(c(1:nrow(completeDNA)), 3, simplify=FALSE)
  pairs <- combn(c(1:nrow(completeDNA)), 2, simplify=FALSE)
  if(method > 1 && length(combos) > 1){
    # Implements the method whereby distance information is used to reject pairs which would likeley be pointless.
    cat("\nTrimming number of triplet comparrisons...")
    if(method==2){                                                    
      binarySeqs <- as.DNAbin(completeDNA)
      distances <- dist.dna(binarySeqs)
      rejectpairs <- pairs[which(distances<0.01)]
      removals <- list()
      for(i in 1:length(combos)){
        for(n in 1:length(rejectpairs)){
          if(all(rejectpairs[[n]] %in% combos[[i]])){
            removals <- append(removals, i)
            break
          }
        }
      }
      combos <- combos[-unlist(removals)]
    } else {
      # Decide pairs to exclude by considering troughs in density distribution.
      if(method==3){
        binarySeqs <- as.DNAbin(completeDNA)
        distances <- dist.dna(binarySeqs)
        distances_density <- density(distances)
        Lows <- cbind(distances_density$x[which(diff(sign(diff(distances_density$y))) == 2)],distances_density$y[which(diff(sign(diff(distances_density$y))) == 2)])
        Lowest <- Lows[which(Lows[,1] == min(Lows[,1])),]
        rejectpairs <- pairs[which(distances < Lowest[1])]
        removals <- list()
        for(i in 1:length(combos)){
          for(n in 1:length(rejectpairs)){
            if(all(rejectpairs[[n]] %in% combos[[i]])){
              removals <- append(removals, i)
              break
            }
          }
        }
        combos <- combos[-unlist(removals)]
      } 
    }
  }
  cat("\nCropping DNA of universally shared sites...")
  # Replaced by much faster impementation: croppedDNA <- completeDNA[,apply(completeDNA,2,function(x) any(c(FALSE,x[-length(x)]!=x[-1])))]
  croppedDNA <- completeDNA[, colSums(completeDNA[-1,] != completeDNA[-nrow(completeDNA), ]) > 0]
  contignames <- rownames(completeDNA)
  cat("\nAll Done!")
  return(as.HybRIDSdna(list(Sequence=completeDNA, CroppedSequence=croppedDNA, Combinations=combos, ContigNames=contignames)))
}