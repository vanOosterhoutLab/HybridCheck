### Functions and classes defining, and for use with sliding window scanning
### functionality in HybridCheck.

## Internal functions.
#' @title Check sliding window settings used for a scan.
#' @description Internal function. Used to make sure the user has entered
#' a sensible option for the sliding window size and step increments.
#' @param winSize A single numeric value. The proposed size of the sliding
#' window in the scan.
#' @param trackLen The A single numeric value. length of the track that the
#' window will travel down.
#' @examples 
#' HybridCheck:::windowSizeChecker(10, 100)
#' HybridCheck:::windowSizeChecker(100, 50)
#' HybridCheck:::windowSizeChecker(100, 10)
#' HybridCheck:::windowSizeChecker(100, 3)
#' @return winSize Single numeric value. The size of the sliding window to
#' be used in the analysis after the check. If an accepable size was chosen by
#' the user, this value will be identical to the input value. 
windowSizeChecker <- function(winSize, trackLen){
  message("\t- Checking the sliding window parameters.")
  if(winSize > trackLen){
    winSize <- as.integer((trackLen / 100) * 10)
    warning("\t\t- The set sliding window size is bigger than the length of the actual informative sites of the contig!")
    message("\t\t- Continuing with analysis but set the sliding window to 10%
            of the sequence length... ")
    message("\t\t- This is equal to ", winSize)
    if(winSize < 1L){
      winSize <- 1L
      message("\t\t- Set the sliding window to 10% of the sequence length, 
              but since this value is below 1, instead setting
              the sliding window length to 1...")
    }
  }
  return(winSize)
}
  

#' @title Generate data-frame with co-ordinates of windows for analysis.
#' @description Internal function. Used to make a data-frame of all the sliding
#' windows: their mid-points, start, ends. This data-frame is used during the 
#' sequence alignment scans and also forms part of the final data table that
#' is produced.
#' @param winSize Single numeric value. The size (in bp) of the sliding window.
#' @param stepSize Single numeric value. The number of base positions the 
#' sliding window jumps along the sequences every time it moves.
#' @param trackLen Single numeric value. Length of the track the window will 
#' travel down.
#' @param bases A numeric vector. The base positions that the sliding window
#' will slide over. Note in many HybridCheck analyses, the sliding window only
#' slides over the informative sites in a sequence alignment i.e. fully
#' conserved, non-informative sites are removed.
#' @return slideFrame A data-frame of six columns. Each column is numeric.
makeWindowFrames <- function(winSize, stepSize, trackLen, bases){
  message("\t- Making all the window frames...")
  if(winSize >= 1L){
    halfWindow <- as.integer(winSize / 2)
    allstepsfrom <- 1 + halfWindow
    allstepsto <- (trackLen - halfWindow) + 1
    allsteps <- seq(from = allstepsfrom, to = allstepsto, by = stepSize)
    windowp1 <- allsteps - halfWindow # All the window start points.
    windowp2 <- allsteps + halfWindow # All the window end points.
    removals <- which(windowp2 > trackLen)
    if(length(removals) > 0) {
      allsteps <- allsteps[-removals]
      windowp1 <- windowp1[-removals]
      windowp2 <- windowp2[-removals]
    }
    slideFrame <- matrix(ncol = 6, nrow = length(windowp1))
    slideFrame[, 1] <- allsteps
    slideFrame[, 2] <- windowp1
    slideFrame[, 3] <- windowp2
    # ActualBP Center
    slideFrame[, 4] <- 
      as.numeric(unlist(lapply(1:length(allsteps),
                               function(i) bases[allsteps[i]])))
    # Actual BP Start
    slideFrame[, 5] <- as.numeric(bases[windowp1])
    # Actual BP End
    slideFrame[, 6] <- as.numeric(bases[windowp2]) 
    # Make placeholders for scan results.
    return(slideFrame)
  } else {
    stop("The sliding window size is less than 1, this is not supposed to be possible.")
  }
}


makeConMats <- function(dnaSequences, pairs){
  conMats <- lapply(pairs, function(i){
    colSums(consensusMatrix(dnaSequences[i]) != 0) > 1
  })
  return(conMats)
}


calculateDistanceTracks <- function(dnaSequences, pairs, resTable){
  # Make consensus matrices for each pair of sequences.
  consensusMatrices <- makeConMats(dnaSequences, pairs)
  mutationCountTracks <- lapply(consensusMatrices, function(x){
    unlist(lapply(seq(nrow(resTable)), function(i){
      stretch <- resTable[i, 2] : resTable[i, 3]
      return(sum(x[stretch]))
    }))
  })
  return(mutationCountTracks)
}

scan.similarity <- function(dna, windowScan, settings){
  message(paste0(" - Scanning sequence similarity for sequences ",
                 paste0(windowScan$SequenceInfo$ContigNames, collapse=", ")))
  cutDNA <- windowScan$SequenceInfo$prepareDNA(dna, ambiguousAreHet)
  windowScan$readAnalysisSettings(settings)
  if(windowScan$SequenceInfo$InformativeUsedLength >= 1){
    windowScan$ScanData$WindowSizeUsed <- 
      windowSizeChecker(windowScan$ScanData$WindowSizeUsed,
                        windowScan$SequenceInfo$InformativeUsedLength)
    if(windowScan$ScanData$WindowSizeUsed >= 1L) {
      Distances <- makeWindowFrames(windowScan$ScanData$WindowSizeUsed,
                                    windowScan$ScanData$StepSizeUsed,
                                    windowScan$SequenceInfo$InformativeUsedLength,
                                    windowScan$SequenceInfo$InformativeUsed)
      message("\t- Scanning Now!")
      tracks <- calculateDistanceTracks(cutDNA, windowScan$SequenceInfo$ContigPairs, Distances)
      Distances <- cbind(Distances, do.call(cbind, tracks))
      colnames(Distances) <- c("WindowCenter", "WindowStart", "WindowEnd", 
                               "ActualCenter", "ActualStart", "ActualEnd",
                               unlist(
                                 lapply(windowScan$SequenceInfo$ContigPairs, 
                                        function(x) paste(x, collapse=":"))))
      dataCols <- seq.int(from = 7, to = ncol(Distances), by = 1)
      Distances[ , dataCols] <- 
        100 - round((as.numeric(Distances[, dataCols]) / 
                       (windowScan$ScanData$WindowSizeUsed + 1)) * 100)
      windowScan$ScanData$Table <- as.data.frame(Distances)
    } else {
      stop("The sliding window size is less than 1, this is not supposed to be possible.")
    }
  } else {
    warning(paste0("There are no informative sites to work on - skipping analysis of sequence combo: ", windowScan$SequenceInfo$ContigNames, collapse = ", "))
  }
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
                FullDNALength = "numeric"
              ),
              
              methods = list(
                initialize =
                  function(seqNames, dnaSequences){
                    if(any(!(seqNames %in% dnaSequences$getSequenceNames()))){
                      stop("Sequence names specified to be scanned are not in 
                           the full set of DNA sequences.")
                    }
                    ContigNames <<- seqNames
                    FullDNALength <<- dnaSequences$getFullLength()
                  },
                
                prepareDNA =
                  function(dna, ambigsAreHet){
                    pulledDNA <- dna$pullDNA(ContigNames)
                    message("\t- Only keeping certain and polymorphic sites.")
                    cutDNA <- cutToInformativeSites(pulledDNA, .self)
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
  

ManyToRefInfo <- 
  setRefClass("ManyToRefInfo",
              
              contains = "SequenceInformation",
              
              fields = list(
                referenceSequence = "character"
                ),
              
              methods = list(
                initialize =
                  function(seqNames, dnaSequence){
                    callSuper(seqNames, dnaSequence)
                    setReferenceSequence(ContigNames[1])
                  },
                
                generateCombinations = 
                  function(){
                    ContigPairs <<- lapply(
                      ContigNames[-which(ContigNames ==
                                           referenceSequence)],
                      function(x){
                        return(c(referenceSequence, x))
                      }
                    )
                  },
                
                setReferenceSequence =
                  function(newReference){
                    if(!(newReference %in% ContigNames)){
                      stop("The chosen reference: ",
                           newReference, " is not in the set of
                                         sequences to be analysed.")
                    }
                    referenceSequence <<- newReference
                    generateCombinations()
                  }
              )
  )


TripletSequenceInfo <- 
  setRefClass("TripletSequenceInfo",
              
              contains = "SequenceInformation",
              
              methods = list(
                initialize = 
                  function(seqNames, dnaSequences){
                    callSuper(seqNames, dnaSequences)
                    generateCombinations()
                  },
                
                generateCombinations =
                  function(){
                    ContigPairs <<- combn(ContigNames, 2, simplify = FALSE)
                  }
              )
  )
    

#' Reference class to store the results from sequence similarity analyses, and
#' block detection runs for a given triplet.
#' @name ScanResults
#' @field TableFile A length 1 character vector storing the temporary filepath
#' of the datatable of sequence similarity table.
#' @field Table Accessor function for the datatable of sequence similarity that
#' is stored on temporary files.
#' @field WindowSizeUsed Single integer value, stores the size of the sliding
#' window used for the scan.
#' @field StepSizeUsed Single integer value, stores the size of the sliding
#' window used for the scan.
ScanResults <- setRefClass("ScanResults",
                              
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
                                    "Method initializes the object, generates
                                    temporary filenames for the sequence 
                                    similarity table."
                                    TableFile <<- tempfile(pattern = "SSTable",
                                                           tmpdir = HCDir)
                                    blankTable()
                                  },
                                
                                blankTable =
                                  function(){
                                    "Method clears the SS analysis results table."
                                    Table <<- data.frame(WindowCenter = NA, 
                                                         WindowStart = NA, 
                                                         WindowEnd = NA,
                                                         ActualCenter = NA, 
                                                         ActualStart = NA, 
                                                         ActualEnd = NA)
                                  },
                                
                                tableIsBlank =
                                  function(){
                                    "Returns TRUE, if the SS analysis table is
                                    blank and no results are contained in it."
                                    return(all(is.na(Table)))
                                  },
                                
                                readSettings =
                                  function(settings){
                                    StepSizeUsed <<- settings$StepSize
                                    WindowSizeUsed <<- settings$WindowSize
                                  },
                                
                                finalize =
                                  function(){
                                    "Called when the object is destroyed, makes
                                    sure to delete the file saved in the system's
                                    temporary directory."
                                    unlink(TableFile)
                                  }
                              )
)
                                
#' Base class used to define specific types of sliding window scans of sequence
#' similarity. Other classes inherit from this class.
#' @name SlidingWindowScan
SlidingWindowScan <- 
  setRefClass("SlidingWindowScan",
              
              fields = list(
                SequenceInfo = "ANY",
                ScanData = "ANY"
              ),
                
              methods = list(
                initialize = 
                  function(HCDir){
                    ScanData <<- ScanResults$new(HCDir)
                  },
                
                noScanPerformed =
                  function(){
                    return(ScanData$tableIsBlank())
                  },
                
                readAnalysisSettings = 
                  function(settingsObj){
                    ScanData$readSettings(settingsObj)
                  },
                
                writeScannedSequences =
                  function(fullDNA){
                    if(noScanPerformed()){
                      stop("Sequences have not been scanned yet.")
                    }
                    fileBaseNames <- paste0(SequenceInfo$ContigNames,
                                            collapse = "_")
                    seqs <- 
                      SequenceInfo$prepareDNAForScan(fullDNA,
                                                     (SequenceInfo$basesResolved() && 
                                                        SequenceInfo$seqsHaveHet()))
                    writeXStringSet(seqs, paste0(fileBaseNames, ".fas"))
                    write(SequenceInfo$InformativeUsed,
                          paste0(fileBaseNames, "_Sites", ".txt", sep = "\n"))
                  },
                
                writeDatedSequences =
                  function(fullDNA){
                    if(noScanPerformed()){
                      stop("Triplet has not been scanned yet.")
                    }
                    fileBaseNames <- paste0(SequenceInfo$ContigNames, collapse = "_")
                    seqs <- 
                      SequenceInfo$prepareDNAForDating(fullDNA,
                                                       (SequenceInfo$basesResolved() &&
                                                          SequenceInfo$seqsHaveHet()))
                    writeXStringSet(seqs, paste0(fileBaseNames, ".fas"))
                  }
              )
  )
                      

ManyToRefScan <- 
  setRefClass("ManyToRefScan",
              
              contains = "SlidingWindowScan",
              
              methods = list(
                initialize =
                  function(sequenceNames, dnaSequences, HCDir){
                    callSuper(HCDir)
                    SequenceInfo <<- ManyToRefInfo$new(sequenceNames,
                                                       dnaSequences)
                  }
              )
  )

TripletScan <- 
  setRefClass("TripletScan",
              
              contains = "SlidingWindowScan",
              
              fields = list(Blocks = "list"),
              
              methods = list(
                initialize =
                  function(sequenceNames, dnaSequences, HCDir){
                    callSuper(HCDir)
                    SequenceInfo <<- TripletSequenceInfo$new(sequenceNames,
                                                             dnaSequences)
                  },
                
                blocksNotFound =
                  function(){
                    "Returns TRUE, if the tables for block are blank and no
                    block findng has been done."
                    return(length(Blocks) == 0)
                  },
                
                blocksNotDated =
                  function(){
                    "Returns TRUE, if the blocks detected have not been tested
                    for significance or had a divergence time estimated."
                    bools <-
                      unlist(lapply(
                        Blocks, 
                        function(y){
                          all(unlist(lapply(y,
                                            function(x){
                                              all(is.na(x$SNPs)) &&
                                                all(is.na(x$CorrectedSNPs)) &&
                                                all(is.na(x$P_Value)) &&
                                                all(is.na(x$P_Threshold)) &&
                                                all(is.na(x$fiveAge)) &&
                                                all(is.na(x$fiftyAge)) &&
                                                all(is.na(x$ninetyFiveAge))
                                            })))
                        }))
                    return(all(bools))
                  },
                
                putativeBlockFind = 
                  function(parameters){
                    "DOCSTRING TO BE COMPLETE"
                    message("\t- Finding blocks for Triplet: ",
                            paste0(SequenceInfo$ContigNames, collapse = ", "))
                    if(noScanPerformed()){
                      stop("No sequence similarity scan data is available for this triplet - can't identify blocks.")
                    }
                    if(parameters$AutoThresholds == TRUE) {
                      message("\t\t- Using the autodetect thresholds method...")
                      message("\t\t- Deciding on suitable thresholds...")
                      thresholds <- autodetect.thresholds(ScanData, parameters)
                    } else {
                      thresholds <- list(parameters$ManualThresholds,
                                         parameters$ManualThresholds, 
                                         parameters$ManualThresholds)
                    }
                    names(thresholds) <-
                      unlist(lapply(
                        combn(SequenceInfo$ContigNames, 2, simplify = F),
                        function(x) paste(x, collapse=":")))
                    message("\t\t- Now beginning Block Search...")
                    Blocks <<- lapply(1:3, function(i){
                      block.find(ScanData$Table[,c(1:6, 6+i)], thresholds[[i]])
                    })
                    names(Blocks) <<- names(thresholds)
                  },
                
                blockDate =
                  function(dnaobj, parameters){
                    "Block Dating method, estimates the ages of blocks detected
                    based on how many mutations are observed in a block and ."
                    message(" - Now dating blocks for sequence triplet ", 
                            paste0(SequenceInfo$ContigNames, sep = ", "))
                    preparedDNA <- SequenceInfo$prepareDNAForDating(dnaobj)
                    ab.blocks <-
                      lapply(Blocks[[1]],
                             function(x){
                               date.blocks(x,
                                           preparedDNA[SequenceInfo$ContigPairs[[1]]],
                                           parameters$MutationRate,
                                           parameters$PValue,
                                           parameters$BonfCorrection,
                                           parameters$DateAnyway,
                                           parameters$MutationCorrection)
                             })
                    ac.blocks <-
                      lapply(Blocks[[2]],
                             function(x){
                               date.blocks(x,
                                           preparedDNA[SequenceInfo$ContigPairs[[2]]],
                                           parameters$MutationRate,
                                           parameters$PValue,
                                           parameters$BonfCorrection,
                                           parameters$DateAnyway,
                                           parameters$MutationCorrection)
                             })
                    bc.blocks <-
                      lapply(Blocks[[3]],
                             function(x){
                               date.blocks(x,
                                           preparedDNA[SequenceInfo$ContigPairs[[3]]],
                                           parameters$MutationRate,
                                           parameters$PValue,
                                           parameters$BonfCorrection,
                                           parameters$DateAnyway,
                                           parameters$MutationCorrection)
                             })
                    out.blocks <- list(ab.blocks, ac.blocks, bc.blocks)
                    names(out.blocks) <- names(Blocks)
                    Blocks <<- out.blocks
                  }
                
              )
  )
              
              
                                                       
                    
                             
                             
                             
                             
                               
                                 
                                   
                                   
                                     
                                                       
                                 
                             











