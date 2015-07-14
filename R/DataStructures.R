#' @title UserBlocks reference class
#' @name UserBlocks
#' @description 
#' The UserBlocks reference class is the class used to store and manipulate user defined blocks in HC.
UserBlocks <- setRefClass("UserBlocks",
                          
                          fields = list(
                            Pairs = "list"
                          ),
                          
                          methods = list(
                            initialize =
                              function(){
                                Pairs <<- list()
                              },
                            
                            hasPairs =
                              function(){
                                return(length(Pairs) > 0)
                              },
                            
                            enforceUserBlocks =
                              function(){
                                if(!hasPairs()){
                                  "Error: UserBlocks object has not been initialized from a HCseq object."
                                }
                              },
                            
                            initializePairsFromDNA =
                              function(dna){
                                if(class(dna) != "HCseq"){stop("Object provided is not of class HCseq")}
                                pairs <- unlist(lapply(combn(unique(dna$getSequenceNames()),2, simplify=F), function(x) paste(x[1], x[2], sep=":")))
                                for (i in pairs){
                                  Pairs[[i]] <<- data.frame(FirstBP=as.numeric(), LastBP=as.numeric(), ApproxBpLength=as.numeric(), SNPs=as.numeric(), CorrectedSNPs=as.numeric(), P_Value=as.numeric(),
                                                            P_Threshold=as.numeric(), fiveAge=as.numeric(), fiftyAge=as.numeric(), ninetyFiveAge=as.numeric())
                                }
                              },
                            
                            addBlock =
                              function(first, last, pair){
                                enforceUserBlocks()
                                index <- processPair(pair)
                                bplength <- abs(last-first)+1
                                Pairs[[index]] <<- rbind(Pairs[[index]], c(first, last, bplength, NA, NA, NA, NA, NA, NA, NA))
                                names(Pairs[[index]]) <<- c("FirstBP", "LastBP", "ApproxBpLength", "SNPs", "CorrectedSNPs", "P_Value", "P_Threshold",
                                                            "fiveAge", "fiftyAge", "ninetyFiveAge")
                              },
                            
                            blankBlocks = function(pair){
                              enforceUserBlocks()
                              index <- processPair(pair)
                              Pairs[[index]] <<- data.frame(FirstBP=as.numeric(), LastBP=as.numeric(), ApproxBpLength=as.numeric(), SNPs=as.numeric(), CorrectedSNPs=as.numeric(), P_Value=as.numeric(),
                                                            P_Threshold=as.numeric(), fiveAge=as.numeric(), fiftyAge=as.numeric(), ninetyFiveAge=as.numeric())
                            },
                            
                            processPair =
                              function(instring){
                                selections <- unlist(strsplit(instring, ":"))
                                if(length(selections) != 2){stop("You must specify two sequences, between which your recombination event occured.")}
                                options <- strsplit(names(Pairs), ":")
                                index <- which(unlist(lapply(lapply(options, function(x) selections %in% x), function(y) all(y))))
                                if(length(index) != 1){stop("Something has gone wrong indexing pairs in triplets - this scenario should not happen, the index of more than or less than one pair should not be possible, contact package maintainer.")}
                                return(index)
                              },
                            
                            dateBlocks =
                              function(sequences, parameters){
                                for(i in 1:length(Pairs)){
                                  if(nrow(Pairs[[i]]) > 0){
                                    pair <- which(sequences$getSequenceNames() %in% unlist(strsplit(names(Pairs)[i], ":")))
                                    Pairs[[i]] <<- date.blocks(Pairs[[i]], sequences, parameters$MutationRate, parameters$PValue, parameters$BonfCorrection, parameters$DateAnyway, parameters$MutationCorrection)
                                  }
                                }
                              },
                            
                            tabulateBlocks = 
                              function(){
                                namesList <- unlist(lapply(1:length(Pairs), function(i) rep(names(Pairs)[[i]], nrow(Pairs[[i]]))))
                                concatTable <- do.call(rbind, Pairs)
                                resultTable <- cbind(namesList, concatTable)
                                colnames(resultTable)[[1]] <- "Sequence_Pair"
                                rownames(resultTable) <- NULL
                                return(resultTable)
                              }
                          ))

#' An internal Reference Class to represent a DNA alignment, read from a FASTA file.
#' @name HCseq
#' @field FullSequence A DNAStringSet containing the full sequence alignment.
#' @field InformativeSequence A DNAStringSet containing the elignment, with uninformative sites removed.
#' @field InformativeBp An integer vector containing the base positions that are informative.
#' @field ReferenceSeq A character vector of length one with the sequence name that is the reference.
#' @field Populations A list of population definitions - a list of vectors containing sequence names.
HCseq <- setRefClass("HCseq",
                     fields = list( 
                       FullSequence = "ANY",
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
                           message("\t- Checking length of sequences...")
                           if(length(unique(width(FullSequence))) > 1){
                             stop("Sequences are not of same length, is this an MSA??")
                           }
                           message(" - Finished DNA input...")
                         },
                       
                       hasDNA =
                         function(){
                           "Returns true if a DNA sequence alignment has been read in and stored in the object. Otherwise returns false."
                           a <- is.initialized(FullSequence)
                           return(a)
                         },
                       
                       enforceDNA =
                         function(){
                           "Enforces some rules about the content of the sequence object and throws errors should they occur."
                           if(!hasDNA()){stop("Error: HCdna object has not got any sequences loaded in.")}
                           #if(length(InformativeSequence) != length(FullSequence)){stop("Error: Number of sequences in the full alignment, and informative alignment are not the same, this shouldn't happen.")}
                         },
                       
                       numberOfSequences =
                         function(){
                           "Returns the number of sequences in the stored alignment."
                           enforceDNA()
                           return(length(FullSequence))
                         },
                       
                       getFullBp =
                         function(){
                           "Returns a vector containing the numbers of the base positions in the aligned sequences."
                           enforceDNA()
                           return(1:getFullLength())
                         },
                       
#                        getInformativeBp =
#                          function(){
#                            "Returns a vector of the base positions of the informative sites in the aligned sequences."
#                            enforceDNA()
#                            return(InformativeBp)
#                          },
                       
                       getFullLength =
                         function(){
                           "Returns the length in base pairs, of the aligned sequences."
                           enforceDNA()
                           return(unique(width(FullSequence)))
                         },
                       
#                        getInformativeLength =
#                          function(){
#                            "Returns the number in base pairs, of informative sites in the aligned sequences."
#                            enforceDNA()
#                            return(unique(width(InformativeSequence)))
#                          },
                       
                       getSequenceNames =
                         function(){
                           "Returns a character vector of the sequence names."
                           enforceDNA()
                           return(names(FullSequence))
                         },
                       
                       pullTriplet =
                         function(selection){
                           "Extracts from the sequence object, a triplet of sequences."
                           enforceDNA()
                           if(length(selection) != 3 || !is.character(selection)){stop("Three sequence names must be provided to pull a triplet of sequences.")}
                           return(FullSequence[selection])
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
                                           #"\nExcluding non-informative sites: ", getInformativeLength(),
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
                                           #" bp</p><p><b>Excluding non-informative sites:</b> ", getInformativeLength(),
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


#' Reference class to store the information from resolving the DNA sequence information in a triplet.
#' @name SequenceInformation
SequenceInformation <- 
  setRefClass("SequenceInformation",
              fields = list(
                ContigNames = "character",
                ContigPairs = "list",
                InformativeUsedLength = "numeric",
                InformativeUsed = "numeric",
                FullDNALength = "numeric",
                NumberOfHet = "numeric",
                Transformations = "data.frame"),
              
              methods = list(
                initialize =
                  function(seqNames, fullLength){
                    ContigNames <<- seqNames
                    ContigPairs <<- combn(ContigNames, 2, simplify = F)
                    FullDNALength <<- fullLength
                    NumberOfHet <<- 0
                    blankTransTable()
                  },
                
                blankTransTable =
                  function(){
                    Transformations <<-
                      data.frame(Base = NA, AmbigOne = NA, AmbigTwo = NA,
                                 AmbigThree = NA, ResolveOne = NA, ResolveTwo = NA,
                                 ResolveThree = NA)
                  },
                
                basesResolved =
                  function(){
                    return(seqsHaveHet() && !all(is.na(Transformations)))
                  },
                
                seqsHaveHet =
                  function(){
                    return(NumberOfHet > 0)
                  },
                
                prepareDNAForScan =
                  function(dna, ambigsAreHet){
                    # Pull the triplet from the DNA structure.
                    seqTriplet <- dna$pullTriplet(ContigNames)
                    stateMatrix <- consensusMatrix(seqTriplet)
                    if(ambigsAreHet){
                      message(" - Treating ambiguous sites of two states as heterozygous.")
                      message(paste0(" - Identifying heterozygous nucleotides for triplet. ",
                                     paste0(ContigNames, collapse = ", ")))
                      ambigTypes <- colSums(as.matrix(stateMatrix[5:10, ]) != 0)
                      NumberOfHet <<- length(which(ambigTypes > 0))
                      if(seqsHaveHet()){
                        if(!basesResolved()){
                          message("   - An appropriate transformation will be decided for each.")
                          message("   - This must be done the first time a triplet is analysed.")
                          heterozygousCode <- list(
                            M = c('A', 'C'), R = c('A', 'G'), W = c('A', 'T'), 
                            S = c('G', 'C'), Y = c('C', 'T'), K = c('G', 'T'))
                          Transformations <<-
                            rbind(transSingleAmb(heterozygousCode, ambigTypes, stateMatrix),
                                  transTwoAmb(heterozygousCode, ambigTypes, stateMatrix),
                                  transThreeAmb(heterozygousCode, ambigTypes, stateMatrix))
                          Transformations <<- Transformations[apply(Transformations, 1,
                                                                    function(x){
                                                                      !all(is.na(x))
                                                                    }), ]
                        }
                        message(" - Triplet of Sequences has heterozygous sites.")
                        message("   - Making alterations to input sequence based on calculated transformations.")
                        seqTriplet <- DNAStringSet(
                          lapply(1:3, function(i) {
                            return(transformSequence(seqTriplet[[i]], Transformations))
                          })
                        ) 
                      }
                    }
                    message("\t- Only keeping certain and polymorphic sites.")
                    conMat <- consensusMatrix(seqTriplet)
                    InformativeUsed <<- 
                      which(
                        (colSums(conMat[c(5:15, 17, 18),]) == 0) &
                          (colSums(conMat != 0) > 1)
                        )
                    InformativeUsedLength <<- length(InformativeUsed)
                    cutDNA <- DNAStringSet(character(length = 3))
                    cutDNA[[1]] <- seqTriplet[[1]][InformativeUsed]
                    cutDNA[[2]] <- seqTriplet[[2]][InformativeUsed]
                    cutDNA[[3]] <- seqTriplet[[3]][InformativeUsed]
                    names(cutDNA) <- names(seqTriplet)
                    return(cutDNA)
                  },
                
                prepareDNAForDating =
                  function(dna, pair){
                    seqTriplet <- dna$pullTriplet(ContigNames)
                    seqNames <- names(seqTriplet)
                    if(seqsHaveHet()){
                      message(" - Making transformations to DNA for dating.")
                      message("  - These transformations are the same ones, used during scan.")
                      seqTriplet <- DNAStringSet(
                        lapply(1:3, function(i){
                          transformSequence(seqTriplet[[i]], Transformations)
                        })
                      )
                    }
                    names(seqTriplet) <- seqNames
                    return(seqTriplet)
                  }
              )
  )


#' Reference class to store the results from sequence similarity analyses, and block detection runs for a given triplet.
#' @name SimilarityScan
#' @field TableFile A length 1 character vector storing the temporary filepath of the datatable of sequence similarity table.
#' @field Table Accessor function for the datatable of sequence similarity that is stored on temporary files.
#' @field WindowSizeUsed Single integer value, stores the size of the sliding window used for the scan.
#' @field StepSizeUsed Single integer value, stores the size of the sliding window used for the scan.
SimilarityScan <- setRefClass("SimilarityScan",
                              
                              fields = list(
                                TableFile = "character",
                                Table = function(value){
                                  if(missing(value)){
                                    read.table(TableFile)
                                  } else {
                                    write.table(value, file = TableFile)
                                  }
                                },
                                WindowSizeUsed = "numeric",
                                StepSizeUsed = "numeric"),
                              
                              methods = list(
                                initialize =
                                  function(HCDir){
                                    "Method initializes the object, generates temporary filenames for the sequence similarity table."
                                    TableFile <<- tempfile(pattern = "SSTable", tmpdir = HCDir)
                                    blankTable()
                                  },
                                
                                blankTable =
                                  function(){
                                    "Method clears the SS analysis results table."
                                    Table <<- data.frame(WindowCenter = NA, WindowStart = NA, WindowEnd = NA,
                                                         ActualCenter = NA, ActualStart = NA, ActualEnd = NA,
                                                         AB = NA, AC = NA, BC = NA)
                                  },
                                
                                tableIsBlank =
                                  function(){
                                    "Returns TRUE, if the SS analysis table is blank and no results are contained in it."
                                    return(all(is.na(Table)))
                                  },
                                
                                finalize =
                                  function(){
                                    "Called when the object is destroyed, makes sure to delete the file saved in the system's temporary directory."
                                    unlink(TableFile)
                                  }
                              )
)

#' Reference class to store and manage triplet data.
#' @name Triplet
#' @field ContigNames character vector of length 3, stores the dna names of which the triplet was made.
#' @field ContigNumbers integer vector of length 3, stores the indecies of the sequences that made the triplet.
Triplet <- setRefClass("Triplet",
                       
                       fields = list(
                         SequenceInfo = "ANY",
                         ContigPairs = "list",
                         ScanData = "ANY",
                         Blocks = "list"
                       ),
                       
                       methods = list(
                         initialize = 
                           function(sequenceNames, fullSeqLength, HCDir){
                             "Initializer function, creates an instance of triplet."
                             SequenceInfo <<- SequenceInformation$new(sequenceNames, fullSeqLength)
                             ScanData <<- SimilarityScan$new(HCDir)
                           },
                         
                         readSettings =
                           function(settings){
                             ScanData$StepSizeUsed <<- settings$StepSize
                             ScanData$WindowSizeUsed <<- settings$WindowSize
                           },
                         
                         noScanPerformed =
                           function(){
                             "Returns TRUE, if the SS analysis table is blank and the informative sites are not known. This is indicative that a scan of the file has not been done yet."
                             return(ScanData$tableIsBlank() && length(SequenceInfo$InformativeUsedLength) == 0)
                           },
                         
                         blocksNotFound =
                           function(){
                             "Returns TRUE, if the tables for block are blank and no block findng has been done."
                             return(length(Blocks) == 0)
                           },
                         
                         blocksNotDated =
                           function(){
                             "Returns TRUE, if the blocks detected have not been tested for significance or had a divergence time estimated."
                             bools <- unlist(lapply(Blocks, function(y) all(unlist(lapply(y, function(x) all(is.na(x$SNPs)) && all(is.na(x$CorrectedSNPs)) && all(is.na(x$P_Value)) && all(is.na(x$P_Threshold)) && all(is.na(x$fiveAge)) && all(is.na(x$fiftyAge)) && all(is.na(x$ninetyFiveAge))   )))))
                             return(all(bools))
                           },
                         
                         # Method for putative block detection.
                         putativeBlockFind = 
                           function(parameters){
                             "DOCSTRING TO BE COMPLETE"
                             message("\t- Finding blocks for Triplet: ", paste0(SequenceInfo$ContigNames, collapse = ", "))
                             if(noScanPerformed()){stop("No sequence similarity scan data is available for this triplet - can't identify blocks.")}
                             if(parameters$AutoThresholds == TRUE) {
                               message("\t\t- Using the autodetect thresholds method...")
                               message("\t\t- Deciding on suitable thresholds...")
                               thresholds <- autodetect.thresholds(ScanData, parameters)
                             } else {
                               thresholds <- list(parameters$ManualThresholds, parameters$ManualThresholds, parameters$ManualThresholds)
                             }
                             names(thresholds) <- unlist(lapply(combn(SequenceInfo$ContigNames, 2, simplify=F), function(x) paste(x, collapse=":")))
                             message("\t\t- Now beginning Block Search...")
                             Blocks <<- lapply(1:3, function(i) block.find(ScanData$Table[,c(1:6, 6+i)], thresholds[[i]]))
                             names(Blocks) <<- names(thresholds)
                           },
                         
                         blockDate =
                           function(dnaobj, parameters){
                             "Block Dating method, estimates the ages of blocks detected based on how many mutations are observed in a block and ."
                             message(" - Now dating blocks for sequence triplet ", 
                                     paste0(SequenceInfo$ContigNames, sep = ", "))
                             preparedDNA <- SequenceInfo$prepareDNAForDating(dnaobj)
                             ab.blocks <- lapply(Blocks[[1]], function(x) date.blocks(x, preparedDNA[SequenceInfo$ContigPairs[[1]]], parameters$MutationRate, parameters$PValue, parameters$BonfCorrection, parameters$DateAnyway, parameters$MutationCorrection))
                             ac.blocks <- lapply(Blocks[[2]], function(x) date.blocks(x, preparedDNA[SequenceInfo$ContigPairs[[2]]], parameters$MutationRate, parameters$PValue, parameters$BonfCorrection, parameters$DateAnyway, parameters$MutationCorrection))
                             bc.blocks <- lapply(Blocks[[3]], function(x) date.blocks(x, preparedDNA[SequenceInfo$ContigPairs[[3]]], parameters$MutationRate, parameters$PValue, parameters$BonfCorrection, parameters$DateAnyway, parameters$MutationCorrection))
                             out.blocks <- list(ab.blocks, ac.blocks, bc.blocks)
                             names(out.blocks) <- names(Blocks)
                             Blocks <<- out.blocks
                           },
                         
                         tabulateBlocks = function(){
                           "Tabulates the blocks for a triplet."
                           message(paste(" - Tabulating blocks for the triplet", 
                                         paste(
                                           SequenceInfo$ContigNames,
                                           collapse=":")))
                           # Check that the tables are present, if they aren't, turn them into blank data.frames.
                           temps <- lapply(1:3, function(i) do.call(rbind, Blocks[[i]]))
                           SS <- lapply(1:3, function(i) floor(as.numeric(rownames(temps[[i]]))))
                           pair <- lapply(1:3, function(i) rep(names(Blocks)[[i]], nrow(temps[[i]])))
                           temp2 <- do.call(rbind, temps)
                           if(blocksNotDated()){
                             temp2 <- cbind(temp2, data.frame(fiveAge = rep(NA, times=nrow(temp2)), fiftyAge = rep(NA, times=nrow(temp2)), ninetyfiveAge = rep(NA, times=nrow(temp2)), SNPnum = rep(NA, times=nrow(temp2)), PValue = rep(NA, times=nrow(temp2)), PThresh = rep(NA, times=nrow(temp2)), MeanAge = rep(NA, times=nrow(temp2)), CorrectedSNPs = rep(NA, times=nrow(temp2))))
                           }
                           temp2["SequencePair"] <- unlist(pair)
                           temp2["SequenceSimilarityThreshold"] <- unlist(SS)
                           if(nrow(temp2) > 0){
                             temp2["Triplet"] <- paste(SequenceInfo$ContigNames, collapse=":")
                           } else {
                             temp2["Triplet"] <- character()
                           }
                           return(temp2)
                         },
                         
                         plotTriplet = function(plottingSettings, begin, end){
                           if("Lines" %in% plottingSettings$What && "Bars" %in% plottingSettings$What){
                             bars <- plotBars(plottingSettings, begin, end)
                             lines <- plotLines(plottingSettings, begin, end)
                             return(arrangeGrob(bars, lines, ncol=1))
                           } else {
                             if("Lines" %in% plottingSettings$What){
                               return(plotLines(plottingSettings, begin, end))
                             }
                             if("Bars" %in% plottingSettings$What){
                               return(plotBars(plottingSettings, begin, end))
                             }
                           }
                         },
                         
                         plotLines =
                           function(plottingSettings, begin, end){
                             "Method plots a lineplot using ggplot2 of the sequence similarity data from the scan."
                             if(noScanPerformed()){stop("No sequence similarity scan has been performed for this triplet.")}
                             combo <- unlist(lapply(combn(SequenceInfo$ContigNames, 2, simplify=FALSE), function(x) paste(x, collapse=":")))
                             data <- ScanData$Table
                             plotting.frame <- data.frame(basepos = rep(data$ActualCenter,3),
                                                          yvalues = c(data$AB, data$AC, data$BC),
                                                          factors = rep(1:3, each = nrow(data)))
                             plot <- ggplot(plotting.frame, aes(x=basepos, y=yvalues)) + geom_line(aes(colour=factor(factors)), show_guide=plottingSettings$Legends, size=0.8) +
                               ylim(0,100) + 
                               xlim(begin, end) +
                               scale_colour_manual(name = "Pairwise Comparrisons", labels=c(combo[1], combo[2], combo[3]),values=c("yellow","purple","cyan")) +
                               xlab("Base Position") +
                               ylab("% Sequence Similarity")
                             plot <- 
                               applyPlottingParams(plot,
                                                   plottingSettings,
                                                   title = paste("Sequence Similarity Between Sequences for Triplet ",
                                                                 SequenceInfo$ContigNames[1],
                                                                 ":",
                                                                 SequenceInfo$ContigNames[2],
                                                                 ":",
                                                                 SequenceInfo$ContigNames[3], sep=""))
                             return(plot)
                           },
                         
                         plotBars =
                           function(plottingSettings, begin, end){
                             "Method plots the heatmap based graphic of bars, from the sequence similarity scan data."
                             if(noScanPerformed()){stop("No sequence similarity scan has been performed for this triplet.")}
                             
                             # Generate the reference colour palette.
                             colourPalette <- expand.grid(A = seq(0, 100, by = 1), B = seq(0, 100, by = 1))
                             colourPalette$RefA <- rgb(green = colourPalette$A, red = 100, blue = colourPalette$B, maxColorValue = 100)
                             colourPalette$RefB <- rgb(green = 100, red = colourPalette$A, blue = colourPalette$B, maxColorValue = 100)
                             colourPalette$RefC <- rgb(green = colourPalette$B, red = colourPalette$A, blue = 100, maxColorValue = 100)
                             
                             # Now figure out the scale and data to go into each vertical bar: TODO - Put this in a function.
                             div <- length(begin:end) / plottingSettings$MosaicScale
                             frame <- data.frame(bpstart = seq(from = begin, to = end, by = div),
                                                 bpend = seq(from = (begin + div) - 1, to = end, by = div))
                             frame$bpX <- round(frame$bpstart +  (div / 2))
                             scanTable <- ScanData$Table
                             AB <- round(apply(frame, 1, function(x) vertbar_create(scanTable, x, 7)))
                             AC <- round(apply(frame, 1, function(x) vertbar_create(scanTable, x, 8)))
                             BC <- round(apply(frame, 1, function(x) vertbar_create(scanTable, x, 9)))
                             rm(scanTable)
                             frame$AB <- AB
                             frame$AC <- AC
                             frame$BC <- BC
                             frame$X <- 1:nrow(frame)
                             rm(AB, AC, BC)
                             if(any(is.nan( as.numeric(frame$AB))) || any(is.nan( as.numeric(frame$AC))) || any(is.nan( as.numeric(frame$BC)))){
                               warning("\nNot a numbers (NaNs)! have been detected in the plotting frame.\n
                                       The common cause of this is a small alignment or few informative sites in the data, 
                                       with a too high MosaicScale parameter.\nThis usually happens at either end of the 
                                       bars and the NaNs will be dealt with my filling them in black.\n\nTo get rid of them use a lower MosaicScale parameter.")
                             }
                             A_mix <- apply(frame, 1, function(x) col_deter(x[c(4,5)], colourPalette[,c(1,2,3)]))
                             B_mix <- apply(frame, 1, function(x) col_deter(x[c(4,6)], colourPalette[,c(1,2,4)]))
                             C_mix <- apply(frame, 1, function(x) col_deter(x[c(5,6)], colourPalette[,c(1,2,5)]))
                             frame$A_mix <- A_mix
                             frame$B_mix <- B_mix
                             frame$C_mix <- C_mix
                             rm(A_mix, B_mix, C_mix)
                             plottingFrame <- data.frame(X = frame$X, Y = rep(c(3, 2, 1), each = plottingSettings$MosaicScale), colour = c(frame$A_mix, frame$B_mix, frame$C_mix))
                             bars <- ggplot(plottingFrame, aes(x = X, y = as.factor(Y))) +
                               geom_raster(aes(fill = colour)) + scale_fill_identity() +
                               xlab("Approximate Base Position") +
                               ylab("Sequence Name") +
                               scale_x_continuous(breaks = c(seq(from = 1, to = plottingSettings$MosaicScale, by = plottingSettings$MosaicScale / 10), plottingSettings$MosaicScale), labels = c(frame$bpX[seq(from = 1, to = plottingSettings$MosaicScale, by = plottingSettings$MosaicScale / 10)], max(frame$bpX))) + 
                               scale_y_discrete(labels = c(SequenceInfo$ContigNames[3], SequenceInfo$ContigNames[2], SequenceInfo$ContigNames[1]))
                             
                             bars <- applyPlottingParams(bars, plottingSettings, title = paste("Sequence Similarity Between Sequences for Triplet ", SequenceInfo$ContigNames[1], ":", SequenceInfo$ContigNames[2], ":", SequenceInfo$ContigNames[3], sep=""))
                             
                             if(plottingSettings$Legends == T){
                               legend <- readPNG(system.file("extdata/rgblegend.png", package="HybridCheck"), TRUE)
                               if (names(dev.cur()) == "windows"){
                                 # windows device doesn’t support semi-transparency so we’ll need
                                 # to flatten the image
                                 transparent <- legend[,,4] == 0
                                 legend <- as.raster(legend[,,1:3])
                                 legend[transparent] <- NA
                               }
                               legendgrob <- grid::rasterGrob(image=legend)
                               bars <- arrangeGrob(bars, legendgrob, widths = c(1, 0.13), ncol = 2)
                             }
                             return(bars)
                             },
                         
                         writeScannedSequences = function(fullDNA){
                           if(noScanPerformed()){stop("Triplet has not been scanned yet.")}
                           fileBaseNames <- paste0(SequenceInfo$ContigNames, collapse = "_")
                           seqs <- SequenceInfo$prepareDNAForScan(fullDNA, (SequenceInfo$basesResolved() && 
                                                                              SequenceInfo$seqsHaveHet()))
                           writeXStringSet(seqs, paste0(fileBaseNames, ".fas"))
                           write(SequenceInfo$InformativeUsed, paste0(fileBaseNames, "_Sites", ".txt", sep = "\n"))
                         },
                         
                         writeDatedSequences = function(fullDNA){
                           if(noScanPerformed()){stop("Triplet has not been scanned yet.")}
                           fileBaseNames <- paste0(SequenceInfo$ContigNames, collapse = "_")
                           seqs <- SequenceInfo$prepareDNAForDating(fullDNA, (SequenceInfo$basesResolved() &&
                                                                                SequenceInfo$seqsHaveHet()))
                           writeXStringSet(seqs, paste0(fileBaseNames, ".fas"))
                         }
                         
                         )
)

vertbar_create <- function(sequenceSimilarityTable, plottingFrameRow, whichComparrison){
  bool1 <- sequenceSimilarityTable$ActualStart <= plottingFrameRow[2]
  bool2 <- plottingFrameRow[1] <= sequenceSimilarityTable$ActualEnd
  index <- which(bool1 == bool2)
  return(mean(sequenceSimilarityTable[index, whichComparrison]))
}

col_deter <- function(invalues, reference){
  if(any(is.nan(as.numeric(invalues)))){
    cols <- "#000000"
  } else {
    cols <- reference[(reference[, 1] == as.numeric(invalues[1])) & reference[,2] == as.numeric(invalues[2]), 3]
  }
  return(cols)
}




#' Reference class storing all triplets in a HC analysis.
#' @name Triplets
#' @description The Triplets reference class stores and manages operations over many Triplet objects.
Triplets <- setRefClass("Triplets",
                        
                        fields = list(
                          triplets = "list"
                        ),
                        
                        methods = list(
                          initialize =
                            function(dna){
                              "Initializes the object."
                              triplets <<- list()
                            },
                          
                          tripletsGenerated =
                            function(){
                              "Returns TRUE if the triplet generations have been prepared for."
                              return(length(triplets) > 0)
                            },
                          
                          matchNames =
                            function(selection){
                              if(!tripletsGenerated()){stop("No triplets have been prepared yet.")}
                              return(unlist(lapply(triplets, function(x) sum(selection %in% x$SequenceInfo$ContigNames))))
                            },
                          
                          deleteAllTriplets =
                            function(){
                              "Removes all current triplet data completely."
                              message(" - Deleting all triplets data.")
                              triplets <<- list()
                            },
                          
                          generateTriplets =
                            function(dna, csettings, basefile){
                              "Initializes all triplet objects based on the combination settings"
                              if(csettings$Modified){
                                if(tripletsGenerated()){
                                  deleteAllTriplets()
                                }
                                message(" - Initializing new triplets data.")
                                seqlength <- dna$getFullLength()
                                seqnames <- dna$getSequenceNames()
                                triplets <<- lapply(csettings$AcceptedCombinations, function(x) Triplet$new(c(x[1], x[2], x[3]), seqlength, basefile))
                                csettings$Modified <- FALSE
                              }
                            },
                          
                          #                           updateTriplets =
                          #                             function(dna, settings, basefile){
                          #                               got <- getAllNames()
                          #                               to.remove <- which(!unlist(lapply(got, function(x) x %in% settings$AcceptedCombinations)))
                          #                               to.add <- which(!unlist(lapply(settings$AcceptedCombinations, function(x) x %in% got)))
                          #                               triplets <<- triplets[-to.remove]
                          #                               triplets <<- append(triplets, makeTriplets(settings$AcceptedCombinations[to.add], dna, basefile))
                          #                             },
                          #                           
                          #                           makeTriplets =
                          #                             function(selection, dna, basefile){
                          #                               seqlength <- dna$getFullLength()
                          #                               return(lapply(selection, function(x) Triplet$new(sequencenumbers = x, sequences = c(x[1], x[2], x[3]), fullseqlength = seqlength, basefile)))
                          #                             },
                          
                          scanTriplets =
                            function(tripletSelections, dna, ambhet, scansettings){
                              if(!tripletsGenerated()){stop("No triplets have been prepared yet.")}
                              tripletsToScan <- getTriplets(tripletSelections)
                              for(tripletToScan in tripletsToScan){
                                scan.similarity(dna, tripletToScan, ambhet, scansettings)
                              }
                            },
                          
                          findBlocks =
                            function(tripletSelections, findSettings){
                              if(!tripletsGenerated()){stop("No triplets have been prepared yet.")}
                              tripletsToFindIn <- getTriplets(tripletSelections)
                              for(triplet in tripletsToFindIn){
                                triplet$putativeBlockFind(findSettings)
                              }
                            },
                          
                          dateBlocks =
                            function(tripletSelections, dateSettings, dna){
                              if(!tripletsGenerated()){stop("No triplets have been prepared yet.")}
                              tripletsToDate <- getTriplets(tripletSelections)
                              for(triplet in tripletsToDate){
                                triplet$blockDate(dna, dateSettings)
                              }
                            },
                          
                          plotTriplets = function(tripletSelections, begin, end, plotSettings){
                            tripletsToPlot <- getTriplets(tripletSelections)
                            return(lapply(tripletsToPlot, function(x) x$plotTriplet(plotSettings, begin, end)))
                          },
                          
                          tabulateBlocks = function(tripletSelections, neat){
                            if(!tripletsGenerated()){stop("No triplets have been prepared yet.")}
                            tripletsToTabulate <- getTriplets(tripletSelections)
                            listedTabulates <- lapply(tripletsToTabulate, function(x) x$tabulateBlocks())
                            collected <- do.call(rbind, listedTabulates)
                            output <- data.frame(collected$Triplet, collected$SequencePair, collected$SequenceSimilarityThreshold, collected$Length,
                                                 collected$First, collected$Last, collected$FirstBP, collected$LastBP, collected$ApproxBpLength, collected$SNPs, collected$CorrectedSNPs, collected$fiveAge, collected$fiftyAge,
                                                 collected$ninetyFiveAge, collected$P_Value, collected$P_Threshold)
                            if(neat){
                              output <- output[,-c(4,5,6)]
                              names(output) <- c("Triplet", "Sequence_Pair","Sequence_Similarity_Threshold","First_BP_Position", "Last_BP_Position", "Approximate_Length_BP", "Number_of_SNPs", "Corrected_Number_of_SNPs", "p=0.05_Age", "p=0.5_Age","p=0.95_Age","P_Value", "P_Thresh")
                            }
                            return(output)
                          },
                          
                          writeTripletsAsScanned = 
                            function(tripletSelections, dna){
                              if(!tripletsGenerated()){stop("No triplets have been prepared yet.")}
                              tripletsToWrite <- getTriplets(tripletSelections)
                              for(triplet in tripletsToWrite){
                                triplet$writeScannedSequences(dna)
                              }
                            },
                          
                          writeTripletsAsDated = 
                            function(tripletsSelections, dna){
                              if(!tripletsGenerated()){stop("No triplets have been prepared yet.")}
                              tripletsToWrite <- getTriplets(tripletSelections)
                              for(triplet in tripletsToWrite){
                                triplet$writeDatedSequences(dna)
                              }
                            },
                          
                          show =
                            function(){
                              cat("A total of ")
                            },
                          
                          tripletCombinations =
                            function(){
                              return(unlist(lapply(triplets, function(x){paste0(x$SequenceInfo$ContigNames, collapse=", ")})))
                            },
                          
                          htmlSummary =
                            function(){
                              "Prints a HTML summary of all the different triplet combinations."
                              tripnames <- paste0(lapply(triplets, function(x){paste0(paste0(x$SequenceInfo$ContigNames, collapse=", "), " <br> ", collapse="")}), collapse="")
                              output <- paste0("<h2>TripletCombinations</h2>")
                              return(tripnames)
                            },
                          
                          getTriplets =
                            function(selections){
                              "Returns a list of references to triplets according to user selection."
                              if(!is.list(selections)){
                                selections <- list(selections)
                              }
                              selections <- unique(selections)
                              if(!is.null(selections) && length(selections) > 0){
                                ind <- numeric()
                                if("ALL" %in% selections){
                                  ind <- 1:length(triplets)
                                } else {
                                  if("NOT.SCANNED" %in% selections){
                                    ind <- c(ind, which(unlist(lapply(triplets, function(x) x$noScanPerformed()))))
                                    selections <- selections[which(selections != "NOT.SCANNED")]
                                  }
                                  if("NOT.SEARCHED" %in% selections){
                                    ind <- c(ind, which(unlist(lapply(triplets, function(x) x$blocksNotFound()))))
                                    selections <- selections[which(selections != "NOT.SEARCHED")]
                                  }
                                  if("NOT.DATED" %in% selections){
                                    ind <- c(ind, which(unlist(lapply(triplets, function(x) x$blocksNotDated()))))
                                    selections <- selections[which(selections != "NOT.DATED")]
                                  }
                                  if(any(unlist(lapply(selections, length)) != 3)){stop("Selections must provide a vector of 3 sequence names.")}
                                  if(any(unlist(lapply(selections, function(x) !is.character(x))))){stop("Selections must be of class character.")}
                                  allNames <- do.call(rbind, getAllNames())
                                  ind <- c(ind, unlist(lapply(selections, function(x) which(allNames[,1] %in% x & allNames[,2] %in% x & allNames[,3] %in% x))))
                                }
                                ind <- unique(ind)
                                return(triplets[ind])
                              } else {
                                stop("No selection of triplet was provided.")
                              }
                            },
                          
                          getAllNames = 
                            function(){
                              "Returns the names of the sequences in each triplet as a list."
                              return(lapply(triplets, function(x) x$SequenceInfo$ContigNames))
                            },
                          
                          getAllIndexes =
                            function(){
                              "Returns the indexes of the sequences (according to their rows in the HC sequence object) in each triplet as a list."
                              return(lapply(triplets, function(x) x$SequenceInfo$ContigIndexes))
                            }
                        )
)


#' A Reference Class for storing and manipulating the results and data from a four taxon test.
#' @name FTTrecord
#' @field P1 Character The population which forms the P1 taxon in the four taxon test.
#' @field P2 Character The population which forms the P2 taxon in the four taxon test.
#' @field P3 Character The population which forms the P3 taxon in the four taxon test.
#' @field A Character The population which forms the ancestral/outgroup taxon in the 
#' four taxon test.
#' 
#' @field numBlocks integer The number of blocks the DNA sequence alignment was split 
#' into in order to perform the test.
#' 
#' @field blockLength integer The number of base pairs to each block the DNA sequence 
#' alignment was split into to perform the test.
#' 
#' @field ABBA numeric The global sum of ABBA sites. 
#' Given by \eqn{(1 - Pr_1) * Pr_2 * Pr_3 * (1 - Pr_4)}. 
#' Where \eqn{Pr_i} is the frequency of the derived allele in the i'th population.
#' 
#' @field ABBA numeric The global sum of BABA sites. 
#' Given by \eqn{Pr_1 * (1 - Pr_2) * Pr_3 * (1 - Pr_4)}.
#' Where \eqn{Pr_i} is the frequency of the derived allele in the i'th population.
#' 
#' @field ABBAcount integer The number of sites for which ABBA was greater than BABA.
#' @field BABAcount integer The number of sites for which BABA was greater than ABBA.
#' 
#' @field globalX2 numeric A chi squared value computed during the four taxon test. 
#' Used to assess whether ABBAcount and BABAcount are significantly different,
#' based on the binomial distribution.
#' 
#' @field X2_P numeric A p-value computer by the Fisher combined probability
#' test based on globalX2. Indicates whether ABBAcount and BABAcount differ significantly.
#' 
#' @field D_jEstimate numeric Patterson's D estimate based on jackknifeing the blocks of data.
#' @field Fd_1DD4_jEstimate numeric A jackknifed estimate of Fd for complete introgression 
#' between populations 2 and 3.
#' 
#' @field Fd_D2D4_jEstimate numeric A jackknifed estimate of Fd for complete introgression 
#' between populations 1 and 3.
#'
#' @field D_jVariance numeric Variance of Pattersons D estimates from jackknife.
#' 
#' @field Fd_1DD4_jVariance numeric Variance of Fd estimates from jackknife.
#' Where Fd is for total introgression between Populations 2 and 3.
#' 
#' @field Fd_D2D4_jVariance numeric Variance of Fd estimates from jackknife.
#' Where Fd is for total introgression between Populations 1 and 3.
#' 
#' @field D_jSD numeric Standard deviation of Patterson's D estimates from jackknife.
#' 
#' @field Fd_1DD4_jSD numeric Standard deviation of Fd estimates from jackknife.
#' Where Fd is for total introgression between Populations 2 and 3.
#' 
#' @field Fd_D2D4_jSD numeric Standard deviation of Fd estimates from jackknife.
#' Where Fd is for total introgression between Populations 1 and 3.
#' 
#' @field D_jZ numeric The Z score of Patterson's D, computed from the jackknife data.
#'  
#' @field Fd_1DD4_jZ numeric The Z score of Fd, computed from the jackknife data.
#' Where Fd is calculated for complete introgression between populations 2 and 3.
#' 
#' @field Fd_D2D4_jZ numeric The Z score of Fd, computed from the jackknife data.
#' Where Fd is calculated for complete introgression between populations 1 and 3.
#' 
#' @field tableFile character A character string indicating the temporary file
#' used to store the dataframe accessed with the table field.
#' 
#' @field table function An accessor function used to access a table of results stored
#' in a temporary file on disk. The location of this file is indicated by the tableFile
#' field.
FTTrecord <- setRefClass("FTTrecord",
                         
                         fields = list(
                           P1 = "character",
                           P2 = "character",
                           P3 = "character",
                           A = "character",
                           numBlocks = "integer",
                           blockLength = "integer",
                           ABBA = "numeric",
                           BABA = "numeric",
                           ABBAcount = "numeric",
                           BABAcount = "numeric",
                           globalX2 = "numeric",
                           X2_P = "numeric",
                           Observed_D = "numeric",
                           Observed_Fd_1DD4 = "numeric",
                           Observed_Fd_D2D4 = "numeric",
                           D_jEstimate = "numeric",
                           Fd_1DD4_jEstimate = "numeric",
                           Fd_D2D4_jEstimate = "numeric",
                           D_jVariance = "numeric",
                           Fd_1DD4_jVariance = "numeric",
                           Fd_D2D4_jVariance = "numeric",
                           D_jSD = "numeric",
                           Fd_1DD4_jSD = "numeric",
                           Fd_D2D4_jSD = "numeric",
                           D_jZ = "numeric",
                           Fd_1DD4_jZ = "numeric",
                           Fd_D2D4_jZ = "numeric",
                           tableFile = "character",
                           table = function(value){
                             if(missing(value)){
                               read.table(tableFile)
                             } else {
                               write.table(value, file = tableFile)
                             }
                           }
                         ),
                         
                         methods = list(
                           initialize =
                             function(p1, p2, p3, a, HCDir){
                               "Initialize the result object."
                               P1 <<- p1
                               P2 <<- p2
                               P3 <<- p3
                               A <<- a
                               tableFile <<- tempfile(pattern = "FTTtable", tmpdir = HCDir)
                               blankTable()
                             },
                           
                           noTestPerformed =
                             function(){
                               "Returns true if a test has not been performed yet."
                               return(all(is.na(table)))
                             },
                           
                           globallySignificant =
                             function(){
                               "Returns true if the test is significant according to the binomial."
                               return((!noTestPerformed()) && (X2_P < 0.05))
                             },
                           
                           blankTable =
                             function(){
                               "Method clears the table."
                               table <<- data.frame(BlockStart = NA, BlockEnd = NA,
                                                    ABBA = NA, BABA = NA, maxABBA_D = NA,
                                                    maxBABA_D = NA, D = NA, Fd = NA)
                             },
                           
                           getPops =
                             function(){
                               "Gets the population names from the result."
                               return(c(P1 = P1, P2 = P2, P3 = P3, A = A))
                             },
                           
                           getTable =
                             function(includeGlobal, neat){
                               "Gets the results table for the blocks used to analyze the data."
                               if(neat && includeGlobal){
                                 out <- cbind(P1, P2, P3, A, table, numBlocks, blockLength,
                                              ABBA, BABA, X2_P, D_jEstimate, Fd_1DD4_jEstimate,
                                              Fd_D2D4_jEstimate, D_jSD, Fd_1DD4_jSD, Fd_D2D4_jSD,
                                              D_jZ, Fd_1DD4_jZ, Fd_D2D4_jZ)
                               } else {
                                 if(includeGlobal){
                                   out <- cbind(P1, P2, P3, A, table, numBlocks, blockLength,
                                                ABBA, BABA, ABBAcount, BABAcount, 
                                                globalX2, X2_P, 
                                                D_jEstimate, Fd_1DD4_jEstimate,
                                                Fd_D2D4_jEstimate, D_jVariance,
                                                Fd_1DD4_jVariance, Fd_D2D4_jVariance,
                                                D_jSD, Fd_1DD4_jSD, Fd_D2D4_jSD,
                                                D_jZ, Fd_1DD4_jZ, Fd_D2D4_jZ)
                                 } else {
                                   out <- cbind(P1, P2, P3, A, table)
                                 }
                               }
                               return(out)
                             }
                         )
)


#' A Reference Class for performing and storing results of four taxon tests for introgression.
#' @name FTTester
#' @field global Logical, whether a global statistic will be calculated for the four taxon tests.
#' @field results A list of reference objects defining the result of a given four taxon test.
FTTester <- setRefClass("FTTester",
                        
                        fields = list(numBlocks = "integer",
                                      taxaCombos = "list",
                                      results = "list"),
                        
                        methods = list(
                          initialize =
                            function(dna){
                              "Initialization method creates object with its default values."
                              numBlocks <<- 50000L
                              results <<- list()
                            },
                          
                          manualTaxaCombos =
                            function(taxas, dna){
                              if(length(taxas) > 0){
                                for(i in taxas){
                                  if(length(unique(i)) != 4){stop("Each taxon combination must provide 4 unique populations: a P1, a P2, a P3 and an A.")}
                                  if(!is.null(names(i))){
                                    inputCheck1 <- all(names(i) %in% c("P1", "P2", "P3", "A"))
                                    if(!inputCheck1){stop("The only names allowed for specifying taxa combos are 'P1', 'P2', 'P3', and 'A'")}
                                  } else {
                                    warning(paste0("No names were provided for the population combination: ", paste0(i, collapse=", "), ".\n",
                                                   "Assuming that P1 is ", i[1], ", that P2 is ", i[2], ", that P3 is ", i[3], " and that A is ", i[4]))
                                    names(i) <- c("P1", "P2", "P3", "A")
                                  }
                                  if(!all(i %in% dna$namesOfPopulations())){stop("You have listed a population name in a taxa combo which does not exist.")}
                                  taxaCombos <<- c(taxaCombos, list(i))
                                }
                              }
                            },
                          
                          autoTaxaCombos =
                            function(dna){
                              if(dna$numberOfPopulations() < 4){
                                stop("Less than 4 populations have been defined - can't form a quartet for an ABBA-BABA test.")
                              }
                              allCombs <- combn(dna$namesOfPopulations(), 4, simplify = FALSE)
                              allDists <- as.matrix(stringDist(dna$FullSequence, method = "hamming")) / dna$getFullLength()
                              generatedCombs <- lapply(allCombs, function(x){
                                out <- list(P1 = NULL, P2 = NULL, P3 = NULL, A = NULL)
                                otus <- x
                                seqsInOtus <- dna$Populations[otus]
                                otuPairs <- combn(otus, 2, simplify = FALSE)
                                distances <- compDist(otuPairs, seqsInOtus, allDists)
                                minOTUs <- distances[which(distances$dist == min(distances$dist)), c("OTU1", "OTU2")]
                                out$P1 <- minOTUs[1,]$OTU1
                                out$P2 <- minOTUs[1,]$OTU2
                                remainingOtus <- otus[which(otus != out$P1 & otus != out$P2)]
                                seqsInOtus2 <- dna$Populations[remainingOtus]
                                P1P2otu <- paste0("(", out$P1, ", ", out$P2, ")")
                                remainingOtus <- c(remainingOtus, P1P2otu)
                                seqsInOtus2[[P1P2otu]] <- unlist(seqsInOtus[c(out$P1, out$P2)])
                                remainingOtuPairs <- combn(remainingOtus, 2, simplify = FALSE)
                                remainingOtuPairs <- remainingOtuPairs[unlist(lapply(remainingOtuPairs, function(x) P1P2otu %in% x))]
                                distances_2 <- compDist(remainingOtuPairs, seqsInOtus2, allDists)
                                minOTUs_2 <- distances_2[which(distances_2$dist == min(distances_2$dist)), c("OTU1", "OTU2")]
                                out$P3 <- minOTUs_2[1, which(minOTUs_2[1,] != P1P2otu)]
                                out$A <- otus[which(!(otus %in% unlist(out)))]
                                checks <- c(
                                  # Check P1 & P2 are closer than P1 and P3.
                                  "d(P1,P2) < d(P1,P3)" = distances[which(apply(distances, 1, function(x){(out$P1 %in% x) & (out$P2 %in% x)})), 3] < 
                                    distances[which(apply(distances, 1, function(x) (out$P1 %in% x) & (out$P3 %in% x))), 3],
                                  # Check P1 & P2 are closer than P1 and A.
                                  "d(P1,P2) < d(P1,PA)" = distances[which(apply(distances, 1, function(x) (out$P1 %in% x) & (out$P2 %in% x))), 3] < 
                                    distances[which(apply(distances, 1, function(x) (out$P1 %in% x) & (out$A %in% x))), 3],
                                  # Check P1 & P2 are closer than P2 and P3.
                                  "d(P1,P2) < d(P2,P3)" = distances[which(apply(distances, 1, function(x){(out$P1 %in% x) & (out$P2 %in% x)})), 3] < 
                                    distances[which(apply(distances, 1, function(x) (out$P2 %in% x) & (out$P3 %in% x))), 3],
                                  # Check P1 & P2 are closer than P2 and A.
                                  "d(P1,P2) < d(P2,A)" = distances[which(apply(distances, 1, function(x) (out$P1 %in% x) & (out$P2 %in% x))), 3] < 
                                    distances[which(apply(distances, 1, function(x) (out$P2 %in% x) & (out$A %in% x))), 3]
                                )
                                if(any(!checks)){
                                  warning(paste0(paste0("Automatically allocating P1, P2, P3 and A designations for population quartet ", paste0(otus, collapse=", "),":\n"),
                                                 paste0("Sanity check: ", names(checks)[which(!checks)], " failed.", collapse = "\n")))
                                }
                                return(out)
                              })
                              manualTaxaCombos(generatedCombs, dna)
                            },
                          
                          hasTaxaCombos =
                            function(){
                              return(length(taxaCombos) > 0)
                            },
                          
                          setNumBlocks =
                            function(value){
                              if(!is.integer(value) || length(value != 1)){stop("Provide only one, integer value.")}
                              numBlocks <<- value
                            },
                          
                          setSettings =
                            function(...){
                              settings <- list(...)
                              parameters <- names(settings)
                              for(i in 1:length(settings)){
                                if(parameters[i] == "numBlocks"){
                                  setNumBlocks(settings[[i]])
                                }
                                if(parameters[i] == "taxaCombos"){
                                  setTaxaCombos(settings[[i]])
                                }
                              }
                            },
                          
                          generateFTTs =
                            function(HCDir){
                              message(" - Initializing new FTtest data.")
                              results <<- lapply(taxaCombos, function(x) FTTrecord$new(x[["P1"]], x[["P2"]], x[["P3"]], x[["A"]], HCDir))
                            },
                          
                          getAllNames = 
                            function(){
                              return(lapply(results, function(x) x$getPops()))
                            },
                          
                          printAllNames =
                            function(){
                              quadNames <- lapply(getAllNames(), function(x){
                                paste0("P1: ", x["P1"], ", P2: ", x["P2"], ", P3: ", x["P3"], ", P4: ", x["A"])
                              })
                              return(quadNames)
                            },
                          
                          getFTTs =
                            function(selections){
                              "Returns a list of references to FTTrecords objects according to user selection."
                              if(!is.list(selections)){
                                selections <- list(selections)
                              }
                              selections <- unique(selections)
                              if(!is.null(selections) && length(selections) > 0){
                                ind <- numeric()
                                if(length(results) != 0){
                                  if("ALL" %in% selections){
                                    ind <- 1:length(results)
                                  } else {
                                    if("NOT.TESTED" %in% selections){
                                      ind <- c(ind, which(unlist(lapply(results, function(x) x$noTestPerformed()))))
                                      selections <- selections[which(selections != "NOT.TESTED")]
                                    }
                                    if("TESTED" %in% selections){
                                      ind <- c(ind, which(unlist(lapply(results, function(x) !x$noTestPerformed()))))
                                      selections <- selections[which(selections != "TESTED")]
                                    }
                                    if("SIGNIFICANT" %in% selections){
                                      ind <- c(ind, which(unlist(lapply(results, function(x) x$globallySignificant()))))
                                      selections <- selections[which(selections != "SIGNIFICANT")]
                                    }
                                    if("PART.SIGNIFICANT" %in% selections){
                                      ind <- c(ind, which(unlist(lapply(results, function(x) x$globallySignificant()))))
                                      selections <- selections[which(selections != "PART.SIGNIFICANT")]
                                    }
                                    if(any(unlist(lapply(selections, length)) != 4)){stop("Selections must provide a vector of 4 sequence names.")}
                                    if(any(unlist(lapply(selections, function(x) !is.character(x))))){stop("Selections must be of class character.")}
                                    allNames <- do.call(rbind, getAllNames())
                                    ind <- c(ind, unlist(lapply(selections, function(x){
                                      which(allNames[,1] %in% x & allNames[,2] %in% x & allNames[,3] %in% x & allNames[,4] %in% x)
                                    })))
                                  }
                                }
                                ind <- unique(ind)
                                return(results[ind])
                              } else {
                                stop("No selection of FTT was provided.")
                              }
                            },
                          
                          runFTTests =
                            function(selections, dna, numBlocks, blocksLen){
                              fttsToTest <- getFTTs(selections)
                              for(ftt in fttsToTest){
                                fourTaxonTest(dna, ftt, numBlocks, blocksLen)
                              }
                            },
                          
                          getResults =
                            function(selections, neat, global){
                              fttsToCollect <- getFTTs(selections)
                              collectedTables <- lapply(fttsToCollect, function(ftt) ftt$getTable(global, neat))
                              return(do.call(rbind, collectedTables))
                            }
                        )
)

