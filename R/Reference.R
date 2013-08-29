# HybRIDS Reference Class Implementation. - Redesign the HybRIDS workflow around one central Reference (R5 Class).

HybRIDS <- setRefClass( "HybRIDS",
                        fields = list( 
                          FullDNAFile = "character",
                          FullDNASequence = function( value ) {
                            if( missing( value ) ){
                              readDNA( FullDNAFile )
                            } else {
                              write.dna( value, file = FullDNAFile, format = "fasta")
                            }
                          },
                          InformativeDNASequence = "matrix",
                          TripletParams = "list",
                          SSAnalysisParams = "list",
                          BlockDetectionParams = "list",
                          BlockDatingParams = "list",
                          Triplets = "list"
                          ),
                        
                         methods = list( initialize = 
                                           function( dnafile = NULL ) {
                                             FullDNAFile <<- tempfile( pattern = "FullDNA" )
                                             TripletParams <<- list(
                                               Method = 1,
                                               SortThreshold = 0.01 )
                                             SSAnalysisParams <<- list(
                                               WindowSize = 100,
                                               StepSize = 1,
                                               SequenceTriplets = list () )
                                             length( BlockDetectionParams ) <<- 4
                                             BlockDetectionParams <<- list( ManualThresholds = c( 90 ), AutoDetectThresholds = TRUE, ManualThresholdFallback = TRUE, SDstringency = 2 )
                                             BlockDatingParams <<- list( MutationRate = 10e-08, RequiredPValue = 0.005 )
                                             if( !is.null( dnafile ) ){
                                               InputDNA( dnafile, format )
                                             }
                                           },
                                         
                                         
                                         InputDNA =
                                           function( intarget, format = "FASTA" ) {
                                             if( format == "FASTA" ){
                                               FullSequence <- InputFasta( intarget )
                                             }
                                             #SequenceNames <<- rownames( CompleteDNA )
                                             #SequenceLength <<- ncol( CompleteDNA )
                                             cat( "\nCropping DNA of universally shared sites..." )
                                             InformativeDNASequence <<- FullSequence[, colSums( FullSequence[-1,] != FullSequence[-nrow( FullSequence ), ] ) > 0]
                                             FullDNASequence <<- FullSequence
                                             #InformativeLength <<- ncol(InformativeDNASequence)
                                             cat( "\nAll Done!" )
                                           }
                                         
                                         
                                         
                                         
#                                          analysis = 
#                                            function() {
#                                              
#                                            }
#                                         # Method for analyzing the sequence similarty of sequences with a sliding window.
#                                         analyzeSS = 
#                                           function() {
#                                             if(length(SequenceCombos) < 2) {
#                                               TripletSubset <- list(CroppedSequence = CroppedDNA, ContigNames = SequenceNames)
#                                               SSTriplets[[1]] <<- seq.similarity(TripletSubset, SSAnalysisParams$WindowSize, SSAnalysisParams$StepSize, SequenceLength, verbose=T)
#                                             } else {
#                                               progress <- txtProgressBar(min=0, max=length(SequenceCombos), style=3)
#                                               for(i in 1:length(SequenceCombos)){
#                                                 setTxtProgressBar(progress, i)
#                                                 Triplet <- SequenceCombos[[i]]
#                                                 TripletSubset <- list( CroppedSequence = CroppedDNA[SequenceCombos[[i]],], ContigNames=SequenceNames[Triplet])
#                                                 SSTriplets[[i]] <<- seq.similarity(TripletSubset, SSAnalysisParams$WindowSize, SSAnalysisParams$StepSize, SequenceLength, verbose=F)
#                                               }
#                                             }
#                                             names(SSTriplets) <<- unlist(lapply(SequenceCombos, function(x) paste(SequenceNames[x[1]], ":", SequenceNames[x[2]], ":", SequenceNames[x[3]], sep="")))
#                                           }
                                          
                          )
                        
                        )

Triplet <- setRefClass( "Triplet",
                          fields = list( DistanceDF = "data.frame",
                                         InformativeDNALength = "numeric",
                                         Contigs = "character",
                                         WindowSizeUsed = "numeric",
                                         StepSizeUsed = "numeric"
                                         ))







                          