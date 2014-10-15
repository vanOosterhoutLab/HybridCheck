#' A Reference Class to represent a DNA alignment, read from a FASTA file.
#'
#' @field FullSequence A Matrix (of Characters) containing the full sequence alignment.
#' @export
HybRIDSseq <- setRefClass("HybRIDSseq",
                            fields = list( 
                              FullSequence = "ANY",
                              InformativeSequence = "ANY",
                              SequenceNames = "character",
                              SequenceLength = "numeric",
                              InformativeLength = "numeric",
                              FullBp = "numeric",
                              InformativeBp = "numeric",
                              NoDNA = "logical"),
                              
                            methods = list( 
                              initialize =
                                function(sequenceInput = NULL, force_Format = NULL) {
                                  NoDNA <<- TRUE
                                  if(!is.null(sequenceInput)){
                                    InputDNA(sequenceInput, force_Format)
                                  }
                                },
                              
                              InputDNA =
                                function(intarget, forceFormat = NULL) {
                                  FullSequence <<- InputSequences(intarget, forceFormat)
                                  FullBp <<- as.numeric(colnames(FullSequence)) # MAKE REDUNDANT
                                  message("Subsetting the informative segregating sites...")
                                  InformativeSequence <<- FullSequence[, sequenceChecker_cpp(FullSequence)] # Cpp code checks for non-informative sites.
                                  InformativeBp <<- as.numeric(colnames(InformativeSequence)) # MAKE REDUNDANT
                                  SequenceLength <<- ncol(FullSequence) # MAKE REDUNDANT
                                  InformativeLength <<- ncol(InformativeSequence) # MAKE REDUNDANT
                                  SequenceNames <<- rownames(FullSequence)
                                  NoDNA <<- FALSE # MAKE REDUNDANT
                                  message("Finished DNA input.")
                                },
                              
                              hasDNA =
                                function(){
                                  a <- is.initialized(FullSequence)
                                  b <- is.initialized(InformativeSequence)
                                  if(a != b){stop("Error: FullSequence is initialized but InformativeSequence is not. This should not happen.")}
                                  return(a)
                                },
                              
                              enforceDNA =
                                function(){
                                  if(!hasDNA()){stop("Error: HybRIDSdna object has not got any sequences loaded in.")}
                                  if(nrow(InformativeSequence) != nrow(FullSequence)){stop("Error: Number of sequences in the full alignment, and informative alignment are not the same, this shouldn't happen.")}
                                },
                              
                              numberOfSequences =
                                function(){
                                  enforceDNA()
                                  return(nrow(InformativeSequence))
                                },
                              
                              getFullBp =
                                function(){
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
                                },
                              
                              plotInf =
                                function(parameters, which="hist") {
                                  if(which != "hist" && which != "bars") stop("You need to specify either 'bars' or 'hist' as the which option!")
                                  if (which == "bars"){
                                    # Figure out the scale and data to go into each vertical bar: TODO - Put this in a function.
                                    div <- FullBp[SequenceLength] / parameters$MosaicScale
                                    bpstart <- seq(from = 1, to = FullBp[SequenceLength], by = div)
                                    bpend <- seq(from=div, to = FullBp[SequenceLength], by = div)
                                    bpX <- round(bpstart + (div / 2))
                                    frame <- data.frame(bpstart = bpstart, bpend = bpend, bpX = bpX)
                                    rm(bpstart, bpend)
                                    Bar <- round(apply(frame, 1, function(x) vertbar_create2(InformativeBp, x)))
                                    plotFrame <- data.frame(bpstart = frame$bpstart, bpend = frame$bpend, X = frame$bpX, numberSNPs = Bar)
                                    thePlot <- ggplot(plotFrame, aes(x = X, y = 1)) + geom_raster(aes(fill = numberSNPs))
                                    thePlot <- applyPlottingParams(thePlot, parameters, title="Number of SNP Polymorphisms Between All Sequences")
                                    thePlot <- thePlot + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
                                  } else {
                                    thePlot <- ggplot(data.frame(InformativeBP = InformativeBp), aes(x = InformativeBP)) + geom_histogram(aes(fill=..count..)) +
                                      scale_fill_gradient2("SNP Count", low="white", high="red")
                                    thePlot <- applyPlottingParams(thePlot, parameters, title="Number of SNP Polymorphisms Between All Sequences")
                                    thePlot <- thePlot + ylab("Count")
                                  }
                                  thePlot <- thePlot + xlab("Base Pair Position")
                                  return(thePlot)
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

sortInput <- function(input, format){
  classOfInput <- class(input)
  if(classOfInput == "character"){
    if(is.null(format)){
      format <- decideFileFormat(input)
    }
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