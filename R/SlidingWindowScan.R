SlidingWindowScan <- setRefClass("SlidingWindowScan",
                       
                       fields = list(
                         SequenceInfo = "ANY",
                         ScanData = "ANY"
                       ),
                       
                       methods = list(
                         initialize = 
                           function(sequenceNames, fullSeqLength, HCDir){
                             "Initializer function, creates an instance of
                             SlidingWindowScan."
                             SequenceInfo <<- 
                               SequenceInformation$new(sequenceNames,
                                                       fullSeqLength)
                             ScanData <<- SimilarityScan$new(HCDir)
                           },
                         
                         noScanPerformed =
                           function(){
                             "Returns TRUE, if the SS analysis table is blank
                             and the informative sites are not known. This is
                             indicative that a scan of the file has not been 
                             done yet."
                             return(ScanData$tableIsBlank() && 
                                      length(SequenceInfo$InformativeUsedLength) == 0)
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