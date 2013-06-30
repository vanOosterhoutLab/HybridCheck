print.HybRIDSdna <- function(x){
  cat("===== HybRIDS DNA sequence object =====\n")
  cat("\nThere are", nrow(x$CroppedSequence), "sequences.\n\n", sep=" ")
  cat("Sequence names:", x$ContigNames, sep="\n")
  cat("\nAligned sequence length:", ncol(x$Sequence), "\n", sep=" ")
  cat("Number of informative sites:", ncol(x$CroppedSequence), "\n", sep=" ")
  cat("Number of triplet analyses to complete:", length(x$Combinations), "\n", sep=" ")
}