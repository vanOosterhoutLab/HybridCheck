# HybRIDS Reference Class Implementation. - Redesign the HybRIDS workflow around one central Reference (R5 Class).

HybRIDS <- setRefClass( "HybRIDS",
                        fields = list( DNAFile = "character",
                                       DNASortMethod = "numeric",
                                       DNASortThreshold = "numeric",
                                       SequenceNames = "character",
                                       SequenceLength = "integer",
                                       CompleteDNA = "matrix",
                                       CroppedDNA = "matrix",
                                       SequenceCombos = "list",
                                       SSAnalysisParams = "list",
                                       SSTriplets = "list",
                                       BlockDetectionParams = "list",
                                       Blocks = "list"
                                       ),
                        
                        methods = list( initialize = 
                                          function( input = NULL, method = 1, RawThresh = 0.01, win.size = 100L, step.size = 1L ) {
                                            DNASortMethod <<- method
                                            DNASortThreshold <<- RawThresh
                                            SSAnalysisParams <<- list( WindowSize = win.size,
                                                                       StepSize = step.size)
                                            BlockDetectionParams <<- list( ManualThresholds = c( 90 ),
                                                                           ThresholdAutodetect = TRUE,
                                                                           FallbackManual = TRUE,
                                                                           SDStringency = 2 )
                                            if( !is.null( input ) ){
                                              DNAFile <<- input
                                              sequenceLoad()
                                            } else {
                                              DNAFile <<- ""
                                            }
                                          },
                                        # Method for analyzing the sequence similarty of sequences with a sliding window.
                                        analyzeSS = 
                                          function() {
                                            if(length(SequenceCombos) < 2) {
                                              TripletSubset <- list(CroppedSequence = CroppedDNA, ContigNames = SequenceNames)
                                              SSTriplets[[1]] <<- seq.similarity(TripletSubset, SSAnalysisParams$WindowSize, SSAnalysisParams$StepSize, SequenceLength, verbose=T)
                                            } else {
                                              progress <- txtProgressBar(min=0, max=length(SequenceCombos), style=3)
                                              for(i in 1:length(SequenceCombos)){
                                                setTxtProgressBar(progress, i)
                                                Triplet <- SequenceCombos[[i]]
                                                TripletSubset <- list( CroppedSequence = CroppedDNA[SequenceCombos[[i]],], ContigNames=SequenceNames[Triplet])
                                                SSTriplets[[i]] <<- seq.similarity(TripletSubset, SSAnalysisParams$WindowSize, SSAnalysisParams$StepSize, SequenceLength, verbose=F)
                                              }
                                            }
                                            names(SSTriplets) <<- unlist(lapply(SequenceCombos, function(x) paste(SequenceNames[x[1]], ":", SequenceNames[x[2]], ":", SequenceNames[x[3]], sep="")))
                                          },
                                        
                                        sequenceLoad =
                                          function() {
                                            dnaSeqLst <- read.fasta(file = DNAFile, 
                                                                    seqtype = "DNA", as.string = FALSE, forceDNAtolower = TRUE,
                                                                    set.attributes = TRUE, legacy.mode = TRUE, seqonly = FALSE, strip.desc = FALSE)
                                            dnaSeqLengths <- sapply(dnaSeqLst, length)
                                            if(all(dnaSeqLengths == dnaSeqLengths[[1]]) == FALSE) stop("\nAll sequences are not of the same length. Please check your input files...\nInput files must be fasta files of an alignment, and 'not' the raw sequence files.\nAborting...")
                                            cat("\nAll sequences are of the same length - good - continuing process...")
                                            cat("\nFormatting data...")
                                            CompleteDNA <<- do.call(rbind, dnaSeqLst)
                                            colnames(CompleteDNA) <<- 1:dnaSeqLengths[[1]]
                                            rownames(CompleteDNA) <<- names(dnaSeqLst)
                                            cat("\nDNA data formatted successfully!")
                                            cat("\nLooking for duplicates...")
                                            dups <- duplicated(CompleteDNA)
                                            if(any(dups)){
                                              cat("\nSome duplicated sequences were found! - We will get rid of these...")
                                              CompleteDNA <<- CompleteDNA[!dups,]
                                            }
                                            if(nrow(CompleteDNA) < 3){
                                              stop("After Removing duplicates, we find that ")
                                            }
                                            SequenceCombos <<- combn(c(1:nrow(CompleteDNA)), 3, simplify=FALSE)
                                            pairs <- combn(c(1:nrow(CompleteDNA)), 2, simplify=FALSE)
                                            if(DNASortMethod > 1 && length(SequenceCombos) > 1){
                                              # Implements the method whereby distance information is used to reject pairs which would likeley be pointless.
                                              cat("\nTrimming number of triplet comparrisons...")
                                              if(DNASortMethod == 2){                                                    
                                                binarySeqs <- as.DNAbin(CompleteDNA)
                                                distances <- dist.dna(binarySeqs)
                                                rejectpairs <- pairs[which(distances < DNASortThreshold)]
                                                removals <- list()
                                                for(i in 1:length(SequenceCombos)){
                                                  for(n in 1:length(rejectpairs)){
                                                    if(all(rejectpairs[[n]] %in% SequenceCombos[[i]])){
                                                      removals <- append(removals, i)
                                                      break
                                                    }
                                                  }
                                                }
                                                SequenceCombos <<- SequenceCombos[-unlist(removals)]
                                              } else {
                                                # Decide pairs to exclude by considering troughs in density distribution.
                                                if(dnaSortMethod == 3){
                                                  binarySeqs <- as.DNAbin(CompleteDNA)
                                                  distances <- dist.dna(binarySeqs, model="raw")
                                                  distances_density <- density(distances)
                                                  Lows <- cbind(distances_density$x[which(diff(sign(diff(distances_density$y))) == 2)],distances_density$y[which(diff(sign(diff(distances_density$y))) == 2)])
                                                  Lowest <- Lows[which(Lows[,1] == min(Lows[,1])),]
                                                  plot(distances_density)
                                                  rejectpairs <- pairs[which(distances < Lowest[1])]
                                                  removals <- list()
                                                  for(i in 1:length(SequenceCombos)){
                                                    for(n in 1:length(rejectpairs)){
                                                      if(all(rejectpairs[[n]] %in% SequenceCombos[[i]])){
                                                        removals <- append(removals, i)
                                                        break
                                                      }
                                                    }
                                                  }
                                                  cat("\nRemoving",length(removals),"triplets")
                                                  SequenceCombos <<- SequenceCombos[-unlist(removals)]
                                                } 
                                              }
                                            }
                                            cat("\nCropping DNA of universally shared sites...")
                                            CroppedDNA <<- CompleteDNA[, colSums(CompleteDNA[-1,] != CompleteDNA[-nrow(CompleteDNA), ]) > 0]
                                            SequenceNames <<- rownames(CompleteDNA)
                                            cat("\nAll Done!")
                                          },
                                        # Method for identifying blocks:
                                        identifyBlocks =
                                          function() {
                                            for(i in 1:length(SSTriplets)){
                                              if(BlockDetectionParams$ThresholdAutodetect == TRUE){
                                                cat("Using the autodetect thresholds method...\n")
                                                cat("Deciding on suitable thresholds...\n")
                                              }
                                              if(BlockDetectionParams$ThresholdAutodetect == TRUE) {
                                                cat("Using the autodetect thresholds method...\n")
                                                cat("Deciding on suitable thresholds...\n")
                                                # Autodetect the thresholds for block identification... uses internal function above.
                                                Thresholds <- autodetect.thresholds(x,SDstringency,manual.thresholds,manualfallback)
                                                names(Thresholds) <- c(paste(x$ContigNames[1],x$ContigNames[2],sep=":"),paste(x$ContigNames[1],x$ContigNames[3],sep=":"),paste(x$ContigNames[2],x$ContigNames[3],sep=":"))
                                              }
                                                
                                                
                                              Blocks[i] <<- identify.blocks(SSTriplets[i], BlockDetectionParams$ManualThresholds, BlockDetectionParams$ThresholdAutodetect, BlockDetectionParams$FallbackManual, BlockDetectionParams$SDStringency)
                                            }
                                            names(Blocks) <<- names(SSTriplets)
                                          }
                                        
                                        
                                            
                                          
                          )
                        
                        )