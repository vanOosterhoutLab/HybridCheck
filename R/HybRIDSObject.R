# HybRIDS Reference Class Implementation. - Redesign the HybRIDS workflow around one central Reference (R5 Class).

HybRIDS <- setRefClass( "HybRIDS",
                        
                        fields = list( 
                          DNA = "HybRIDSseq",
                          TripletParams = "list",
                          SSAnalysisParams = "list",
                          BlockDetectionParams = "list",
                          BlockDatingParams = "list",
                          LastTripletSelection = "numeric",
                          PlottingParams = "list",
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
                                             PlottingParams <<- list( What = c("Bars", "Lines"), Title = TRUE, CombinedTitle = FALSE, 
                                                                      TitleSize = 14, XTicks = TRUE, YTicks = TRUE, 
                                                                      XLabel = TRUE, YLabel = TRUE, Legends = TRUE )
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
                                               Triplets <<- lapply( SSAnalysisParams$TripletCombinations, function(x) HybRIDStriplet$new( sequences = c( DNA$SequenceNames[x[1]], DNA$SequenceNames[x[2]], DNA$SequenceNames[x[3]] ), fullseqlength = DNA$SequenceLength ) )
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
                                           function( Selections = "all" ) {
                                             if( !is.character( Selections ) ) stop( "option 'which' must be 'all' or a vector of the sequence triplets you want to use e.g. 'Seq1:Seq2:Seq3'" )
                                             if( Selections == "all" ) {
                                               if( length( SSAnalysisParams$TripletCombinations ) < 2 ) {
                                                 cat( "Only one triplet to analyze the sequence similarity of..." )
                                                 seq.similarity( DNA$InformativeSequence, Triplets[[1]], SSAnalysisParams$WindowSize, SSAnalysisParams$StepSize, DNA$SequenceLength, DNA$InformativeBp, verbose = T )
                                               } else {
                                                 cat( "Analyzing the sequence similarity of all the triplets...\n" )
                                                 progress <- txtProgressBar( min = 0, max = length(SSAnalysisParams$TripletCombinations), style = 3 )
                                                 for( i in 1:length( SSAnalysisParams$TripletCombinations ) ) {
                                                   setTxtProgressBar( progress, i )
                                                   seq.similarity( DNA$InformativeSequence[ SSAnalysisParams$TripletCombinations[[i]], ], Triplets[[i]], SSAnalysisParams$WindowSize, SSAnalysisParams$StepSize, DNA$SequenceLength, DNA$InformativeBp, verbose = F )
                                                 }
                                               }
                                             } else {
                                               indexTriplets( Selections )
                                               for( i in LastTripletSelection ){
                                                 cat( "Now analysing sequence similarity of triplet", unlist(SSAnalysisParams$TripletCombinations[i]), "\n" )
                                                 seq.similarity( DNA$InformativeSequence[ unlist(SSAnalysisParams$TripletCombinations[[i]]), ], Triplets[[i]], SSAnalysisParams$WindowSize, SSAnalysisParams$StepSize, DNA$SequenceLength, DNA$InformativeBp, verbose = F )
                                               } 
                                             }
                                           },
                                         
                                         # GGplot method for HybRIDS object - activates submethods of triplets.
                                         plotSS =
                                           function( Selections, What = c( "Lines", "Bars" ), Combine = TRUE, ... ) {
                                             if( !is.character( Selections ) ) stop( "option 'Selections' must be a vector of the sequence triplets you want to use e.g. 'Seq1:Seq2:Seq3'" )
                                             Parameters <- list( ... )
                                             for( i in Selections ) {
                                               cat("Selection", i)
                                               if( length( unlist( strsplit(i, ":") ) ) == 3 ) {
                                                 indexTriplets( i )
                                                 if( "Lines" %in% What && !"Bars" %in% What ) {
                                                   outplot <- Triplets[[LastTripletSelection]]$plotLines( ... )
                                                 }
                                                 if( !"Lines" %in% What && "Bars" %in% What ) {
                                                   outplot <- Triplets[[LastTripletSelection]]$plotBars( ... )
                                                 }
                                                 if( "Lines" %in% What && "Bars" %in% What ) {
                                                   outplot <- arrangeGrob( Triplets[[LastTripletSelection]]$plotBars( ... ),
                                                                Triplets[[LastTripletSelection]]$plotLines( ... ),
                                                                ncol = 1 )
                                                 }
                                                 return(outplot)
                                               } else {
                                                 if( length( unlist( strsplit( i, ":" ) ) ) == 2 && Combine == TRUE ) {
                                                   indexTriplets( i )
                                                   if("LabelFontSize" %in% names(Parameters)){
                                                     if(!is.numeric(Parameters$LabelFontSize)) stop("Parameter LabelFontSize must be a number.")
                                                     LabelFontSize <- Parameters$LabelFontSize
                                                   } else {
                                                     LabelFontSize <- 12
                                                   }
                                                   if("LegendFontSize" %in% names(Parameters)){
                                                     if(!is.numeric(Parameters$LegendFontSize)) stop("Parameter LegendFontSize must be a number.")
                                                     LegendFontSize <- Parameters$LabelFontSize
                                                   } else {
                                                     LegendFontSize <- 12
                                                   }
                                                   if( "Title" %in% names( Parameters ) ) {
                                                     if( !is.logical( Parameters$Title ) ) stop( "Parameter Title must be logical." )
                                                     Title <- Parameters$Title
                                                   } else {
                                                     Title <- TRUE
                                                   }
                                                   if( Title == TRUE && ("TitleSize" %in% names( Parameters ) || "CombinedTitle" %in% names(Parameters))) {
                                                     if(!is.numeric(Parameters$TitleSize) || length(Parameters$TitleSize) > 1) stop("Parameter TitleSize must be a numeric values")
                                                     TitleSize <- Parameters$TitleSize
                                                   } else {
                                                     TitleSize <- 12
                                                   }
                                                   if("CombinedTitle" %in% names(Parameters)){
                                                     if(!is.logical(Parameters$CombinedTitle)) stop("Parameter CombinedTitle must be logical.")
                                                     CombinedTitle <- Parameters$CombinedTitle
                                                   } else {
                                                     CombinedTitle <- FALSE
                                                   }
                                                   if("TickSize" %in% names(Parameters)){
                                                     if(!is.numeric(Parameters$TickSize)) stop("Parameter TickSize must be a number.")
                                                     TickSize <- Parameters$TickSize
                                                   } else {
                                                     TickSize <- 12
                                                   }
                                                   if("TickColour" %in% names(Parameters)){
                                                     if(!is.character(Parameters$TickColour)) stop("Parameter TickSize must be a character string e.g. \"red\".")
                                                     TickColour <- Parameters$TickColour
                                                   } else {
                                                     TickColour <- "black"
                                                   }
                                                   if( "Lines" %in% What ){
                                                     dflength <- sum( unlist( lapply( Triplets[LastTripletSelection], function(x) nrow(x$SSTable) ) ) )
                                                     plotting.frame <- data.frame( matrix( nrow = dflength, ncol = 9 ) )
                                                     names(plotting.frame) <- c("WindowCenter", "WindowStart", "WindowEnd", "ActualCenter", "ActualStart", "ActualEnd", "SSVals", "TripletSet", "xvals")
                                                     plotting.frame$xvals <- unlist( lapply( Triplets[LastTripletSelection], function(x) 1:nrow( x$SSTable ) ) )
                                                     plotting.frame$WindowCenter <- unlist( lapply( Triplets[LastTripletSelection], function(x) x$SSTable$WindowCenter ) )
                                                     plotting.frame$WindowStart <- unlist( lapply( Triplets[LastTripletSelection], function(x) x$SSTable$WindowStart ) )
                                                     plotting.frame$WindowEnd <- unlist( lapply( Triplets[LastTripletSelection], function(x) x$SSTable$WindowEnd ) )
                                                     plotting.frame$ActualCenter <- unlist( lapply( Triplets[LastTripletSelection], function(x) x$SSTable$ActualCenter ) )
                                                     plotting.frame$ActualStart <- unlist( lapply( Triplets[LastTripletSelection], function(x) x$SSTable$ActualStart ) )
                                                     plotting.frame$ActualEnd <- unlist( lapply( Triplets[LastTripletSelection], function(x) x$SSTable$ActualEnd ) )
                                                     plotting.frame$SSVals <- unlist( lapply( Triplets[LastTripletSelection], function(x) x$returnPair( unlist( strsplit( Selections, ":" ) )[1], unlist( strsplit( Selections, ":" ) )[2] ) ) )
                                                     plotting.frame$TripletSet <- as.factor( unlist( lapply( Triplets[LastTripletSelection], function(x) rep( paste( c(x$ContigNames[1], x$ContigNames[2], x$ContigNames[3]), collapse=":" ), nrow( x$SSTable ) ) ) ) )
                                                     outplotLines <- ggplot( plotting.frame, aes( x = ActualCenter, y = SSVals ) ) +
                                                       geom_line( aes( colour = TripletSet ), show_guide = T, size = 0.8 ) +
                                                       ylab( "% Sequence Similarity" ) +
                                                       xlab( "Base Position" ) +
                                                       theme( 
                                                         title = element_text(size = 14, colour = "black", face = "bold" ),
                                                         axis.title.y=element_text(size=LabelFontSize),
                                                         axis.text.y=element_text(size=TickSize, colour=TickColour),
                                                         axis.text.x=element_text(size=TickSize, colour=TickColour),
                                                         axis.title.x=element_text(size=LabelFontSize))
                                                     if(Title == TRUE){
                                                       outplotLines <- outplotLines + ggtitle(paste("Sequence similarity for sequence pair ", i, " in all triplets in which it occurs", sep="" ))
                                                     }
                                                   }
                                                   if( "Bars" %in% What ){
                                                     if( "Mosaic.Scale" %in% names(Parameters) ) {
                                                       if(!is.integer(Parameters$Mosaic.Scale) || length(Parameters$Mosaic.Scale)) stop("The Mosaic.Scale Parameter must be an integer (e.g. 500L)")
                                                       Mosaic.Scale <- Parameters$Mosaic.Scale
                                                     } else {
                                                       Mosaic.Scale <- 500L
                                                     }
                                                     bars <- lapply( Triplets[LastTripletSelection], function(x) x$plotBars( exportDat = T, ... ) )
                                                     pairs <- unlist( lapply( Triplets[LastTripletSelection], function(x) x$returnPair( unlist( strsplit( i, ":" ) )[1], unlist( strsplit( i, ":" ) )[2], data = F ) ) )
                                                     datasize <- sum( unlist( lapply( bars, function(x) nrow(x) ) ) )
                                                     plotting.frame2 <- data.frame( matrix( nrow = datasize, ncol = 3 ) )
                                                     names( plotting.frame2 ) <- c("X","Y","SequenceSimilarity")
                                                     plotting.frame2$SequenceSimilarity <- unlist( lapply( 1:length( bars ), function(i) if( pairs[i] == 1 ){
                                                       bars[[i]]$AB
                                                     } else {
                                                       if( pairs[i] == 2 ){
                                                         bars[[i]]$AC
                                                       } else {
                                                         if( pairs[i] == 3 ){
                                                           bars[[i]]$BC
                                                         }
                                                       }
                                                     } ) )
                                                     plotting.frame2$X <- unlist( lapply( bars, function(x) x$X ) )
                                                     contignames <- unlist( lapply( Triplets[LastTripletSelection], function(x) paste( x$ContigNames[1], ":", x$ContigNames[2], ":", x$ContigNames[3], sep = "") ) )
                                                     plotting.frame2$Y <- rep(1:length(bars), times = unlist(lapply(bars, function(x) nrow(x))))
                                                     bpX <- bars[[1]]$bpX
                                                     yaxislab <- unlist(lapply(Triplets[LastTripletSelection], function(x) paste(x$ContigNames[1],x$ContigNames[2],x$ContigNames[3],sep=":")))
                                                     outplotBars <- ggplot( plotting.frame2, aes( x = X, y = as.factor(Y) ) ) +
                                                       geom_raster( aes( fill = SequenceSimilarity ) ) +
                                                       xlab( "Base Position" ) +
                                                       scale_x_continuous( breaks = c(seq( from = 1, to = Mosaic.Scale, by = Mosaic.Scale / 10 ), Mosaic.Scale), labels = c(bpX[seq( from = 1, to = Mosaic.Scale, by = Mosaic.Scale / 10 )], max(bpX)) ) +
                                                       scale_y_discrete( labels = as.character(yaxislab) ) +
                                                       scale_fill_gradient2(high="red", low="blue", midpoint=33.3) +
                                                       theme( 
                                                         title = element_text(size = TitleSize, colour = "black", face = "bold" ),
                                                         axis.title.y=element_blank(),
                                                         axis.text.y=element_text(size=TickSize, colour=TickColour),
                                                         axis.text.x=element_text(size=TickSize, colour=TickColour),
                                                         axis.title.x=element_text(size=LabelFontSize))
                                                     if( Title == TRUE || CombinedTitle == TRUE ){
                                                       outplotBars <- outplotBars + ggtitle(paste("Sequence similarity for sequence pair ", i, " in all triplets in which it occurs", sep="" ))
                                                     }
                                                   }
                                                   if( "Lines" %in% What && !"Bars" %in% What ) {
                                                     return(outplotLines)
                                                   } else {
                                                     if( !"Lines" %in% What && "Bars" %in% What ) {
                                                       return(outplotBars)
                                                     } else {
                                                       if( "Lines" %in% What && "Bars" %in% What) {
                                                         return( arrangeGrob( outplotBars, outplotLines, ncol = 1 ) )
                                                       }
                                                     }
                                                   }
                                                 }
                                               }
                                             }
                                             },
                                           
                                         
                                         # Method for indexing triplets.
                                         indexTriplets =
                                           function( selections ) {
                                             if( !is.character( selections ) ) stop( "option 'which' must be a vector of the sequence triplets you want to use e.g. 'Seq1:Seq2:Seq3'" )
                                             processedSelections <- strsplit( selections, split=":" )
                                             threes <- processedSelections[which( lapply( processedSelections, function(x) length(x) ) == 3 )]
                                             twos <- processedSelections[which( lapply( processedSelections, function(x) length(x) ) == 2 )]
                                             ones <- processedSelections[which( lapply( processedSelections, function(x) length(x) ) == 1 )]
                                             threes <- lapply( threes, function(x) which( DNA$SequenceNames %in% x ) )
                                             twos <- lapply( twos, function(x) which( DNA$SequenceNames %in% x ) )
                                             ones <- lapply( ones, function(x) which( DNA$SequenceNames %in% x ) )
                                             threes <- which( SSAnalysisParams$TripletCombinations %in% threes )
                                             if( length(twos) > 0 ) {
                                               twos <- which( unlist( lapply( twos, function(y) lapply( SSAnalysisParams$TripletCombinations, function(x) all( y %in% x ) ) ) ) )
                                             } else {
                                               twos <- c()
                                             }
                                             if( length(ones) > 0 ) {
                                               ones <- which( unlist( lapply( ones, function(y) lapply( SSAnalysisParams$TripletCombinations, function(x) any( y %in% x ) ) ) ) )
                                             } else {
                                               ones <- c()
                                             }
                                             # Now let's get rid of redunancies and assign the selection to LastTripletSelection.
                                             LastTripletSelection <<- unique( c( threes, twos, ones ) )
                                           }
                                         



                                          
                          )
                        
                        )











                          