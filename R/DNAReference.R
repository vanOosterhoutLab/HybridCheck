#' A Reference Class to represent a DNA alignment, read from a FASTA file.
#' @name HybRIDSseq
#' @import methods
#' @field FullSequence A Matrix (of Characters) containing the full sequence alignment.
#' @field InformativeSequence A matrix (of Characters) containing the elignment, with uninformative sites removed.
HybRIDSseq <- setRefClass("HybRIDSseq",
                            fields = list( 
                              FullSequence = "ANY",
                              InformativeSequence = "ANY"),
                              
                            methods = list( 
                              initialize =
                                function(sequenceInput = NULL) {
                                  "Initializes the object, may be provided with a filepath to a sequence file, currently only FASTA is supported."
                                  if(!is.null(sequenceInput)){
                                    InputDNA(sequenceInput)
                                  }
                                },
                              
                              InputDNA =
                                function(intarget) {
                                  "Reads in sequences from file and appropriately modifies fields of the object."
                                  FullSequence <<- InputSequences(intarget)
                                  message("Subsetting the informative segregating sites...")
                                  InformativeSequence <<- FullSequence[, sequenceChecker_cpp(FullSequence)] # Cpp code checks for non-informative sites.
                                  message("Finished DNA input.")
                                },
                              
                              hasDNA =
                                function(){
                                  "Returns true if a DNA sequence alignment has been read in and stored in the object. Otherwise returns false."
                                  a <- is.initialized(FullSequence)
                                  b <- is.initialized(InformativeSequence)
                                  if(a != b){stop("Error: FullSequence is initialized but InformativeSequence is not. This should not happen.")}
                                  return(a)
                                },
                              
                              enforceDNA =
                                function(){
                                  "Enforces some rules about the content of the sequence object and throws errors should they occur."
                                  if(!hasDNA()){stop("Error: HybRIDSdna object has not got any sequences loaded in.")}
                                  if(nrow(InformativeSequence) != nrow(FullSequence)){stop("Error: Number of sequences in the full alignment, and informative alignment are not the same, this shouldn't happen.")}
                                },
                              
                              numberOfSequences =
                                function(){
                                  "Returns the number of sequences in the stored alignment."
                                  enforceDNA()
                                  return(nrow(InformativeSequence))
                                },
                              
                              getFullBp =
                                function(){
                                  "Returns the "
                                  enforceDNA()
                                  return(as.numeric(colnames(FullSequence)))
                                },
                              
                              getInformativeBp =
                                function(){
                                  enforceDNA()
                                  return(as.numeric(colnames(InformativeSequence)))
                                },
                              
                              getFullLength =
                                function(){
                                  enforceDNA()
                                  return(ncol(FullSequence))
                                },
                              
                              getInformativeLength =
                                function(){
                                  enforceDNA()
                                  return(ncol(InformativeSequence))
                                },
                              
                              getSequenceNames =
                                function(){
                                  enforceDNA()
                                  return(rownames(InformativeSequence))
                                },
                              
                              pullTriplet =
                                function(selection){
                                  if(length(selection) != 3 || !is.character(selection)){stop("Three sequence names must be provided to pull a triplet of sequences.")}
                                  return(InformativeSequence[selection, ])
                                }
                            ))


# INTERNAL FUNCTIONS:

InputSequences <- function(input, Format) {
  dna <- sortInput(input, Format)
  dna <- checkForDuplicates(dna)
  dna <- as.character(dna)
  colnames(dna) <- 1:ncol(dna)
  message("Done...")
  return(dna)
}

decideFileFormat <- function(input){
  if(grepl(".fas", input)){
    message("File to be read is expected to be FASTA format...")
    return("fasta")
  }
}

sortInput <- function(input){
  classOfInput <- class(input)
  if(classOfInput == "character"){
    format <- decideFileFormat(input)
    message("Reading in sequence file...")
    dna <- read.dna(file=input, format=format, as.matrix=TRUE)
  } else {
    if(classOfInput == "DNAbin"){
      message("Class of input is DNAbin from ape package.")
      dna <- input
    } else {
      error("Input is not a valid path to a DNA file, nor is it a valid DNA object, for example, DNAbin from package ape.")
    }
  }
  return(dna)
}

checkForDuplicates <- function(dna){
  message("Looking for duplicates (sequences with p_distances of 0)...")
  distances <- dist.dna(dna, model="N")
  if(any(distances == 0)){
    message("Duplicated sequences were found! - These will be deleted...")
    indicies <- distrowcol(which(distances == 0), attr(distances, "Size"))
    dna <- dna[-indicies[,2],]
    message("Double Checking for duplicated sequences again to be safe...")
    distances <- dist.dna(dna, model="N")
    if(any(distances == 0)) stop("Duplicates were still found in the sequences - this should not happen - aborting. Inform package maintainer.")
  }
  return(dna)
}

distrowcol <- function(ix,n){
  nr <- ceiling(n-(1+sqrt(1+4*(n^2-n-2*ix)))/2)
  nc <- n-(2*n-nr+1)*nr/2+ix+nr
  return(cbind(nr,nc))
}

is.initialized <- function(x){
  return(class(x) != "uninitializedField")
}