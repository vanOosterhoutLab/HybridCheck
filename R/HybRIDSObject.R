# HybRIDS Reference Class Implementation. - Redesign the HybRIDS workflow around one central Reference (R5 Class).

HybRIDS <- setRefClass( "HybRIDS",
                        
                        fields = list( 
                          DNA = "HybRIDSseq",
                          TripletParams = "list",
                          SSAnalysisParams = "list",
                          BlockDetectionParams = "list",
                          BlockDatingParams = "list",
                          Triplets = "list"
                          ),
                        
                         methods = list( initialize = 
                                           function( dnafile = NULL, format = "FASTA" ) {
                                             TripletParams <<- list(
                                               Method = 1,
                                               SortThreshold = 0.01 )
                                             SSAnalysisParams <<- list(
                                               WindowSize = 100,
                                               StepSize = 1,
                                               TripletCombinations = list () )
                                             length( BlockDetectionParams ) <<- 4
                                             BlockDetectionParams <<- list( ManualThresholds = c( 90 ), AutoThresholds = TRUE, ManualFallback = TRUE, SDstringency = 2 )
                                             BlockDatingParams <<- list( MutationRate = 10e-08, PValue = 0.005 )
                                             if( !is.null( dnafile ) ){
                                               DNA$InputDNA( dnafile, format )
                                             }
                                           },
                                         
                                         # Method for generating the Triplet Combinations...
                                         makeTripletCombos =
                                           function() {
                                             if( nrow( DNA$InformativeSequence ) < 3 ){
                                               stop( "After Removing duplicated there are not enough sequences to make any triplets...\nWe can't do much with these sequences...\nAborting.\n" )
                                             } else {
                                               SSAnalysisParams$TripletCombinations <<- combn( c( 1:nrow( DNA$InformativeSequence ) ), 3, simplify=FALSE )
                                               pairs <- combn( c( 1:nrow( DNA$InformativeSequence ) ), 2, simplify=FALSE )
                                               if( TripletParams$Method > 1 && length( combos ) > 1 ) {
                                                 # Implements the method whereby distance information is used to reject pairs which would likeley be pointless.
                                                 cat( "\nTrimming number of triplet comparrisons..." )
                                                 if( TripletParams$Method == 2 ) {                                                    
                                                   distances <- dist.dna( as.DNAbin( DNA$FullSequence ), model = "raw" )
                                                   rejectpairs <- pairs[ which( distances < RawThresh ) ]
                                                   removals <- list()
                                                   for( i in 1:length( TripletCombinations ) ) {
                                                     for( n in 1:length( rejectpairs ) ) {
                                                       if( all( rejectpairs[[n]] %in% TripletCombinations[[i]] ) ) {
                                                         removals <- append( removals, i )
                                                         break
                                                       }
                                                     }
                                                   }
                                                   SSAnalysisParams$TripletCombinations <<- SSAnalysisParams$TripletCombinations[ -unlist( removals ) ]
                                                 } else {
                                                   # Decide pairs to exclude by considering troughs in density distribution.
                                                   if( TripletParams$Method == 3 ){
                                                     distances <- dist.dna( as.DNAbin( DNA$FullSequence ), model="raw" )
                                                     distances_density <- density( distances )
                                                     Lows <- cbind( distances_density$x[ which( diff( sign( diff( distances_density$y ) ) ) == 2 ) ], distances_density$y[ which( diff( sign( diff( distances_density$y ) ) ) == 2 ) ] )
                                                     Lowest <- Lows[ which( Lows[ ,1 ] == min( Lows[ ,1 ] ) ) , ]
                                                     plot( distances_density )
                                                     rejectpairs <- pairs[ which( distances < Lowest[ 1 ] ) ]
                                                     removals <- list()
                                                     for( i in 1:length( TripletCombinations ) ) {
                                                       for( n in 1:length( rejectpairs ) ) {
                                                         if( all( rejectpairs[[n]] %in% TripletCombinations[[i]] ) ) {
                                                           removals <- append( removals, i )
                                                           break
                                                         }
                                                       }
                                                     }
                                                     cat( "\nRemoving",length( removals ),"triplets" )
                                                     SSAnalysisParams$TripletCombinations <<- SSAnalysisParams$TripletCombinations[ -unlist( removals ) ]
                                                   }
                                                 }
                                               }
                                               Triplets <<- lapply( SSAnalysisParams$TripletCombinations, function(x) HybRIDStriplet$new( sequences = c( DNA$SequenceNames[x[1]], DNA$SequenceNames[x[2]], DNA$SequenceNames[x[3]] ) ) )
                                               names( Triplets ) <<- unlist( lapply( SSAnalysisParams$TripletCombinations, function(x) paste( DNA$SequenceNames[x[1]], DNA$SequenceNames[x[2]], DNA$SequenceNames[x[3]], sep = ":" ) ) )
                                             }
                                           },
                                         
                                         # Method for setting any parameter for any stage.
                                         setParameters =
                                           function( Step, ... ) {
                                             if( Step != "SSAnalysis" && Step != "BlockDetection" && Step != "BlockDating" ){
                                               stop("You need to specify a valid analysis 'Step' to alter the paramerters of.\nThe steps are SSAnalysis, BlockDetection, and BlockDating.")
                                             }
                                             Parameters <- list( ... )
                                             if( Step == "SSAnalysis" ) {
                                               if( "WindowSize" %in% names( Parameters ) ) {
                                                 if( !is.numeric( Parameters$WindowSize ) ) stop( "The Window Size parameter for Sequence Similarity analysis must be numeric..." )
                                                 if( length( Parameters$WindowSize) > 1 ) stop( " The Window Size parameter for Sequence Similarity analysis must be a single number..." )
                                                 SSAnalysisParams$WindowSize <<- Parameters$WindowSize
                                               }
                                               if( "StepSize" %in% names( Parameters ) ) {
                                                 if( !is.numeric( Parameters$StepSize ) ) stop( "The Step Size parameter for Sequence Similarity analysis must be numeric..." )
                                                 if( length( Parameters$StepSize) > 1 ) stop( " The Step Size parameter for Sequence Similarity analysis must be a single number..." )
                                                 SSAnalysisParams$StepSize <<- Parameters$StepSize
                                               }
                                             }
                                             if( Step == "BlockDetection" ) {
                                               if( "ManualThresholds" %in% names( Parameters ) ) {
                                                 if( !is.numeric( Parameters$ManualThresholds ) ) stop( "The Manual Threshold parameter for block detection must be a vector of numeric values..." )
                                                 BlockDetectionParams$ManualThresholds <<- Parameters$ManualThresholds
                                               }
                                               if( "AutoThresholds" %in% names( Parameters ) ) {
                                                 if( !is.logical( Parameters$AutoThresholds ) ) stop( "The AutoThresholds parameter for block detection must be logical (TRUE or FALSE)..." )
                                                 if( length( Parameters$AutoThresholds ) > 1 ) stop( "The AutoThresholds parameter for block detection must be a single logical value (TRUE or FALSE)..." )
                                                 BlockDetectionParams$AutoThresholds <<- Parameters$AutoThresholds
                                               }
                                               if( "ManualFallBack" %in% names( Parameters ) ) {
                                                 if( !is.logical( Parameters$ManualFallBack ) ) stop( "The ManualFallBack parameter for block detection must be logical (TRUE or FALSE)..." )
                                                 if( length( Parameters$ManualFallBack ) > 1 ) stop( "The ManualFallBack parameter for block detection must be a single logical value (TRUE or FALSE)..." )
                                                 BlockDetectionParams$ManualFallBack <<- Parameters$ManualFallBack
                                               }
                                               if( "SDstringency" %in% names( Parameters ) ) {
                                                 if( !is.numeric( Parameters$SDstringency ) ) stop( "The SDstringency parameter for block detection must be a numeric value..." )
                                                 if( length( Parameters$SDstringency ) > 1 ) stop( "The SDstringency parameter for block detection must be a single numeric value..." )
                                                 BlockDetectionParams$ManualThresholds <<- Parameters$ManualThresholds
                                               }
                                             }
                                             if( Step == "BlockDating" ) {
                                               if( "MutationRate" %in% names( Parameters ) ) {
                                                 if( !is.numeric( Parameters$MutationRate ) ) stop( "The Mutation Rate parameter for block dating must be a numeric value..." )
                                                 if( length( Parameters$MutationRate ) > 1 ) stop( "The Mutation Rate parameter for block dating must be a single numeric value..." )
                                                 BlockDatingParams$MutationRate <<- Parameters$MutationRate
                                               }
                                               if( "PValue" %in% names( Parameters ) ) {
                                                 if( !is.numeric( Parameters$PValue ) ) stop( "The P Value threshold parameter for block dating must be a numeric value..." )
                                                 if( length( Parameters$PValue ) > 1 ) stop( "The P Value threshold parameter for block dating must be a single numeric value..." )
                                                 BlockDatingParams$PValue <<- Parameters$PValue
                                               }
                                             }
                                           },
                                         
                                         # Method for analyzing the sequence similarity of triplets of sequences.
                                         analyzeSS = 
                                           function( triplet = "all" ) {
                                             if( !is.character( triplet ) ) stop( "option 'which' must be 'all' or a vector of the sequence triplets you want to use e.g. 'Seq1:Seq2:Seq3'" )
                                             if( triplet == "all" ) {
                                               if( length( SSAnalysisParams$TripletCombinations ) < 2 ) {
                                                 cat( "Only one triplet to analyze\nAnalyzing triplet of sequences" )
                                                 seq.similarity( DNA$InformativeSequence, Triplets[[1]], SSAnalysisParams$WindowSize, SSAnalysisParams$StepSize, DNA$SequenceLength, DNA$InformativeBp, verbose = T )
                                               } else {
                                                 progress <- txtProgressBar( min = 0, max = length( SSAnalysisParams$TripletCombinations ), style = 3 )
                                                 for( i in 1:length( SSAnalysisParams$TripletCombinations ) ) {
                                                   setTxtProgressBar( progress, i )
                                                   seq.similarity( DNA$InformativeSequence[ SSAnalysisParams$TripletCombinations[[i]], ], Triplets[[i]], SSAnalysisParams$WindowSize, SSAnalysisParams$StepSize, DNA$SequenceLength, DNA$InformativeBp, verbose = F )
                                                 }
                                               }
                                             } else {
                                               for( i in triplet ) {
                                                 index <- which( unlist( lapply( SSAnalysisParams$TripletCombinations, function(x) all( which( DNA$SequenceNames %in% unlist( strsplit( i, split=":" ) ) ) %in% x ) ) ) )
                                                 cat( "Analyzing Triplet:", SSAnalysisParams$TripletCombinations[[index]] )
                                                 seq.similarity( DNA$InformativeSequence[ SSAnalysisParams$TripletCombinations[[index]], ], Triplets[[index]], SSAnalysisParams$WindowSize, SSAnalysisParams$StepSize, DNA$SequenceLength, DNA$InformativeBp, verbose = F )
                                               }
                                             }
                                           }



                                          
                          )
                        
                        )











                          