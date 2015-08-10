### Manipulation of DNA sequences.

## General functions that may be useful outside of their internal use in
## HybridCheck.

# Return alignment of only the selected base positions.
extractBases <- function(dna, indexes){
  if(!is(dna, "DNAStringSet")){stop("Argument dna needs to be of class DNAStringSet")}
  if(!is.numeric(indexes)){stop("Argument indexes needs to be of class Numeric")}
  i <- IntegerList(rep.int(list(indexes), length(dna)))
  newSeq <- dna[i]
  return(newSeq)
}


# Find inforative sites in DNAStringSet and return a DNAStringSet
# containing only the informative sites.
cutToInformativeSites <- function(sequences, sequenceInfo){
  conMat <- consensusMatrix(sequences)
  sequenceInfo$InformativeUsed <- 
    which(
      (colSums(conMat[c(5:15, 17, 18),]) == 0) &
        (colSums(conMat != 0) > 1)
    )
  sequenceInfo$InformativeUsedLength <- 
    length(sequenceInfo$InformativeUsed)
  return(extractBases(sequences, sequenceInfo$InformativeUsed))
}


#' An internal Reference Class to represent a DNA alignment, read from a FASTA file.
#' @name HCseq
#' @field FullSequence A DNAStringSet containing the full sequence alignment.
#' @field InformativeSequence A DNAStringSet containing the elignment, with 
#' uninformative sites removed.
#' @field InformativeBp An integer vector containing the base positions that 
#' are informative.
#' @field ReferenceSeq A character vector of length one with the sequence name 
#' that is the reference.
#' @field Populations A list of population definitions - a list of vectors 
#' containing sequence names.
#' @field referenceSeq Single string. The name of the sequence to use as 
#' the reference, when doing sliding window sequence similarity scans of many 
#' sequences against one reference.
#' @field manyToRefTable An instance of the SimilarityScan Reference class
#' stores the data-frame of results from sliding window sequence similarity
#' scans of many sequences against one reference, to a temporary file.
#' @field winSizeUsed A single numberic value representing the size (in base
#' pairs)of the sliding window used in sliding window sequence similarity 
#' scans of many sequences against one reference.
#' @field stepSizeUsed A single numberic value representing the size of the 
#' sliding window used in sliding window sequence similarity scans of many 
#' sequences against one reference.
HCseq <- setRefClass("HCseq",
                     fields = list( 
                       FullSequence = "ANY",
                       Populations = "list",
                       ScanAnalysis = "ANY"
                     ),
                     
                     methods = list( 
                       initialize =
                         function(sequenceInput = NULL, hybCheckDir){
                           "Initializes the object, may be provided with a
                           filepath to a sequence file, currently only FASTA
                           is supported."
                           
                           if(!is.null(sequenceInput)){
                             InputDNA(sequenceInput, hybCheckDir)
                           }
                         },
                       
                       InputDNA =
                         function(intarget, hybCheckDir){
                           "Reads in sequences from file and appropriately
                           modifies fields of the object."
                           FullSequence <<- sortInput(intarget)
                           FullSequence <<- checkForDuplicates(FullSequence)
                           message("\t- Checking length of sequences...")
                           if(length(unique(width(FullSequence))) > 1){
                             stop("Sequences are not of same length, is this an MSA??")
                           }
                           message(" - Finished DNA input...")
                           ScanAnalysis <<-
                             ManyToRefScan$new(getSequenceNames(),
                                               .self,
                                               hybCheckDir)
                         },
                       
                       hasDNA =
                         function(){
                           "Returns true if a DNA sequence alignment has been
                           read in and stored in the object. Otherwise returns false."
                           a <- is.initialized(FullSequence)
                           return(a)
                         },
                       
                       enforceDNA =
                         function(){
                           "Enforces some rules about the content of the 
                           sequence object and throws errors should they occur."
                           if(!hasDNA()){
                             stop("Error: HCdna object has not got any sequences loaded in.")
                           }
                         },
                       
                       numberOfSequences =
                         function(){
                           "Returns the number of sequences in the stored
                           alignment."
                           enforceDNA()
                           return(length(FullSequence))
                         },
                       
                       getFullBp =
                         function(){
                           "Returns a vector containing the numbers of the base
                           positions in the aligned sequences."
                           enforceDNA()
                           return(1:getFullLength())
                         },
                       
                       getFullLength =
                         function(){
                           "Returns the length in base pairs, of the aligned 
                           sequences."
                           enforceDNA()
                           return(unique(width(FullSequence)))
                         },
                       
                       getSequenceNames =
                         function(){
                           "Returns a character vector of the sequence names."
                           enforceDNA()
                           return(names(FullSequence))
                         },
                       
                       pullDNA =
                         function(selection){
                           "Extracts from the sequence object, a triplet of
                           sequences."
                           enforceDNA()
                           return(FullSequence[selection])
                         },
                       
                       setPopulations =
                         function(pops){
                           "Define which sequences form a population. Provide a
                           list of vectors containing either integers or
                           sequence names."
                           enforceDNA()
                           if(length(pops) > 0){
                             pops <- popIntegersToNames(pops, getSequenceNames())
                           }
                           Populations <<- pops
                           if(is.null(names(Populations))){
                             names(Populations) <<- paste("unnamed", 1:length(Populations), sep = "_")
                           } else {
                             for(i in 1:length(Populations)){
                               ind <- 1
                               if(names(Populations)[[i]] == ""){
                                 names(Populations)[[i]] <<- paste("unnamed", ind, sep = "_")
                                 ind <- ind + 1
                               }
                             }
                           }
                         },
                       
                       oneSeqOnePop = function(){
                         "Function assigns one population per sequence."
                         enforceDNA()
                         setPopulations(as.list(getSequenceNames()))
                       },
                       
                       numberOfPopulations =
                         function(){
                           "Returns the number of populations assigned."
                           return(length(Populations))
                         },
                       
                       hasPopulations =
                         function(){
                           "Returns TRUE when a set of populations has been
                           defined. Otherwise returns FALSE."
                           return(length(Populations) >= 1)
                         },
                       
                       namesOfPopulations =
                         function(){
                           "Returns the names of the populations."
                           return(names(Populations))
                         },
                       
                       textSummary =
                         function(){
                           "Creates a character vector of the summary of the
                           sequence object."
                           start <- paste0("DNA Sequence Information:\n",
                                           "-------------------------\nAn alignment of ", numberOfSequences(), 
                                           " sequences.\n\nFull length of alignment: ", getFullLength(),
                                           "\n\nSequence names:\n")
                           names <- getSequenceNames()
                           end <- paste0(lapply(1:length(names), function(i) paste0(i, ": ", names[i])), collapse = "\n")
                           if(length(Populations) == 0){
                             pops <- "No populations have been specified."
                           } else {
                             pops <- paste0(lapply(1:length(Populations), 
                                                   function(i){paste0(namesOfPopulations()[i],
                                                                      ": ",
                                                                      paste0(Populations[[i]], collapse = ", ")
                                                   )}), collapse = "\n")
                           }
                           return(paste0(start, end, "\n\nPopulations:\n", pops, "\n"))
                         },
                       
                       htmlSummary =
                         function(){
                           "Creates a character vector of the summary of the
                           sequence object, formatted as HTML."
                           start <- paste0("<h2>DNA Sequence Information:</h2>",
                                           "<p>An alignment of ", 
                                           numberOfSequences(),
                                           " sequences.</p><p><b>Full length of alignment:</b> ",
                                           getFullLength(),
                                           " bp</p><p><b>Sequence names:</b><br>")
                           names <- getSequenceNames()
                           end <- paste0(lapply(1:length(names), function(i) paste0(i, ": ", names[i])), collapse="<br>")
                           if(length(Populations) == 0){
                             pops <- "No populations have been specified."
                           } else {
                             pops <- paste0(lapply(1:length(Populations),
                                                   function(i){paste0(namesOfPopulations()[i], ": ",
                                                                      paste0(Populations[[i]], 
                                                                             collapse = ", "))}),
                                            collapse = "<br>")
                           }
                           return(paste0(start, end, "<br><br><b>Populations:</b><br>", pops, "<br>"))
                         },
                       
                       show =
                         function(){
                           "Prints a text summary of the object to console."
                           cat(textSummary())
                         },
                       
                       scanManyToReference =
                         function(reference = NULL, winSize = 100L, stepSize = 1L){
                           "Performs a sliding window scan of sequence
                           similarity which compares many sequences to a
                           single reference."
                           if(!is.null(reference)){
                             ScanAnalysis$SequenceInfo$setReference(reference)
                           }
                           settings <- SSAnalysisSettings$new(winSize = winSize, stepSize = stepSize)
                           scan.similarity(.self, ScanAnalysis, settings)
                         }
                       )
                     )

## Internal functions specifically for the HCseq reference class. 

sortInput <- function(input){
  classOfInput <- class(input)
  message(" - Reading in sequence file...")
  if(classOfInput == "character"){
    dna <- readDNAStringSet(filepath = input, format = "fasta")
  } else {
    if(classOfInput == "DNAStringSet"){
      message("\t- Class of input is DNAStringSet from Biostrings package...")
      dna <- input
    } else {
      if(classOfInput == "DNAbin"){
        message("\t- Class of input is DNAbin from the ape package...")
        dna <- DNAStringSet(unlist(apply(as.character(input), 1, function(x) paste0(x, collapse = ""))))
      } else {
        stop("Input is not a valid path to a DNA file, nor is it a valid DNA object.")
      }
    }
  }
  return(dna)
}

checkForDuplicates <- function(dna){
  message("\t- Looking for duplicates...")
  duplicates <- duplicated(dna)
  if(any(duplicates)){
    message("\t- Duplicated sequences were found! - These will be deleted...")
    dna <- dna[which(!duplicates)]
  }
  return(dna)
}

is.initialized <- function(x){
  return(class(x) != "uninitializedField")
}





