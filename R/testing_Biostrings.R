library(Biostrings)
origMAlign <- readDNAStringSet(filepath = "~/Dropbox/MySequences.fas", format = "fasta")
# Make a consensus matrix which counts the number of occurences of the bases at that location.
consensusM <- consensusMatrix(origMAlign)
polymorphic <- which(colSums(consensusM != 0) > 1)
newalign <- DNAStringSet()
length(newalign) <- 1
for(i in 1:length(origMAlign))
  newalign[[i]] <- origMAlign[[i]][polymorphic]

