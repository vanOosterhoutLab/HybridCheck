#' An internal Reference Class to represent DNA data, read from a FASTA file.
#' @name SequenceData
#' @field FullSequence A DNAStringSet containing the full sequence alignment.
#' @field InformativeSequence A DNAStringSet containing the elignment, with uninformative sites removed.
#' @field InformativeBp An integer vector containing the base positions that are informative.
#' @field ReferenceSeq A character vector of length one with the sequence name that is the reference.
#' @field Populations A list of population definitions - a list of vectors containing sequence names.
SequenceData <- setRefClass("SequenceData",
                            fields = list(
                              FullSequence = "ANY",
                              InformativeSequence = "ANY",
                              InformativeBp = "integer",
                              ReferenceSeq = "character",
                              Populations = "list"
                              ),

                            methods = list(
                              initialize =
                                function(sequenceInput = NULL){
                                  "Initializes the object, may be provided with a filepath to a sequence file, currently only FASTA is supported."
                                  if(!is.null(sequenceInput)){
                                    InputDNA(sequenceInput)
                                  }
                                },

                              InputDNA =
                                function(intarget){
                                  "Reads in sequences from file and appropriately modifies fields of the object."
                                  FullSequence <<- sortInput(intarget)
                                  FullSequence <<- checkForDuplicates(FullSequence)
                                  message("Checking length of sequences...")
                                  if(length(unique(width(FullSequence))) > 1){
                                    stop("Sequences are not of same length, is this an MSA??")
                                  }
                                  message("Subsetting the informative segregating sites...")
                                  consensusM <- consensusMatrix(FullSequence)
                                  notUnknown <- consensusM[15, ] == 0
                                  polymorphic <- colSums(consensusM != 0) > 1
                                  InformativeBp <<- which(notUnknown & polymorphic)
                                  index <- rep.int(list(InformativeBp), length(FullSequence))
                                  InformativeSequence <<- FullSequence[index]
                                  names(InformativeSequence) <<- names(FullSequence)
                                  ReferenceSeq <<- names(FullSequence)[1]
                                  message("Finished DNA input...")
                                },

                              addBAMSAMs =
                                function(files){

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
                                  if(!hasDNA()){stop("Error: HCdna object has not got any sequences loaded in.")}
                                  if(length(InformativeSequence) != length(FullSequence)){stop("Error: Number of sequences in the full alignment, and informative alignment are not the same, this shouldn't happen.")}
                                },

                              numberOfSequences =
                                function(){
                                  "Returns the number of sequences in the stored alignment."
                                  enforceDNA()
                                  return(length(InformativeSequence))
                                },

                              getFullBp =
                                function(){
                                  "Returns a vector containing the numbers of the base positions in the aligned sequences."
                                  enforceDNA()
                                  return(1:getFullLength())
                                },

                              getInformativeBp =
                                function(){
                                  "Returns a vector of the base positions of the informative sites in the aligned sequences."
                                  enforceDNA()
                                  return(InformativeBp)
                                },

                              getFullLength =
                                function(){
                                  "Returns the length in base pairs, of the aligned sequences."
                                  enforceDNA()
                                  return(unique(width(FullSequence)))
                                },

                              getInformativeLength =
                                function(){
                                  "Returns the number in base pairs, of informative sites in the aligned sequences."
                                  enforceDNA()
                                  return(unique(width(InformativeSequence)))
                                },

                              getSequenceNames =
                                function(){
                                  "Returns a character vector of the sequence names."
                                  enforceDNA()
                                  return(names(InformativeSequence))
                                },

                              pullTriplet =
                                function(selection){
                                  "Extracts from the sequence object, a triplet of sequences."
                                  enforceDNA()
                                  if(length(selection) != 3 || !is.character(selection)){stop("Three sequence names must be provided to pull a triplet of sequences.")}
                                  return(InformativeSequence[selection])
                                },

                              setPopulations =
                                function(pops){
                                  "Define which sequences form a population. Provide a list of vectors containing either integers or sequence names."
                                  enforceDNA()
                                  if(length(pops) > 0){
                                    if(any(!unlist(lapply(pops, function(x) is.integer(x) || is.character(x))))){stop("Need to provide a list of groups of sequence names or integers representing sequence numbers.")}
                                    pops <- lapply(pops, function(x){
                                      if(is.integer(x)){
                                        return(getSequenceNames()[x])
                                      } else {
                                        return(x)
                                      }
                                    })
                                    if(any(table(unlist(pops)) > 1)){stop("Entered a sequence name or number in more than one group.")}
                                    if(any(!unlist(lapply(pops, function(x) all(x %in% getSequenceNames()))))){stop("Some sequences specified in the populations are not in the sequence data.")}
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
                                  "Returns TRUE when a set of populations has been defined. Otherwise returns FALSE."
                                  return(length(Populations) >= 1)
                                },

                              namesOfPopulations =
                                function(){
                                  "Returns the names of the populations."
                                  return(names(Populations))
                                },

                              textSummary =
                                function(){
                                  "Creates a character vector of the summary of the sequence object."
                                  start <- paste0("DNA Sequence Information:\n",
                                                  "-------------------------\nAn alignment of ", numberOfSequences(),
                                                  " sequences.\n\nFull length of alignment: ", getFullLength(),
                                                  "\nExcluding non-informative sites: ", getInformativeLength(),
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
                                  "Creates a character vector of the summary of the sequence object, formatted as HTML."
                                  start <- paste0("<h2>DNA Sequence Information:</h2>",
                                                  "<p>An alignment of ", numberOfSequences(),
                                                  " sequences.</p><p><b>Full length of alignment:</b> ", getFullLength(),
                                                  " bp</p><p><b>Excluding non-informative sites:</b> ", getInformativeLength(),
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
                                }
                            ))


# INTERNAL FUNCTIONS:

sortInput <- function(input){
  classOfInput <- class(input)
  message("Reading in sequence file...")
  if(classOfInput == "character"){
    dna <- readDNAStringSet(filepath = input, format = "fasta")
  } else {
    if(classOfInput == "DNAStringSet"){
      message("Class of input is DNAStringSet from Biostrings package...")
      dna <- input
    } else {
      if(classOfInput == "DNAbin"){
        message("Class of input is DNAbin from the ape package...")
        dna <- DNAStringSet(unlist(apply(as.character(input), 1, function(x) paste0(x, collapse = ""))))
      } else {
        stop("Input is not a valid path to a DNA file, nor is it a valid DNA object.")
      }
    }
  }
  return(dna)
}

checkForDuplicates <- function(dna){
  message("Looking for duplicates...")
  duplicates <- duplicated(dna)
  if(any(duplicates)){
    message("Duplicated sequences were found! - These will be deleted...")
    dna <- dna[which(!duplicates)]
  }
  return(dna)
}

is.initialized <- function(x){
  return(class(x) != "uninitializedField")
}
