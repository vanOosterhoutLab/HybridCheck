#' @docType package
#' ...
#' @import ape ggplot2 grid gridExtra png


#' HybRIDS reference class
#' @export
HybRIDS <- setRefClass( "HybRIDS",
                        
                        fields = list( 
                          DNA = "ANY",
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
                                             PlottingParams <<- list( What = c("Bars", "Lines"), PlotTitle = TRUE, CombinedTitle = FALSE, 
                                                                      TitleSize = 14, TitleFace="bold", TitleColour = "black", XLabels = TRUE, YLabels = TRUE,
                                                                      XTitle = TRUE, XTitleFontSize = 12, XTitleColour = "black", XLabelSize = 10, XLabelColour = "black",
                                                                      YTitle = TRUE, YTitleFontSize = 12, YTitleColour = "black", YLabelSize = 10, YLabelColour = "black",
                                                                      Legends = TRUE, LegendFontSize = 12, MosaicScale = 500)
                                             DNA <<- HybRIDSseq$new()
                                             if( !is.null( dnafile ) ){
                                               DNA$InputDNA( dnafile, format="FASTA")
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
                                         
                                         # Method for displaying parameters.
                                         showParameters = 
                                           function( Step = "ALL") {
                                             if( Step!= "TripletGeneration" && Step != "SSAnalysis" && Step != "BlockDetection" && Step != "BlockDating" && Step != "Plotting" && Step != "ALL"){
                                               stop("You need to specify a valid analysis 'Step' to alter the paramerters of.\nThe steps are TripletGeneration, SSAnalysis, BlockDetection, BlockDating, and Plotting.\n")
                                             }
                                             if( Step == "SSAnalysis"){
                                               cat("Parameters for the Sliding Window Sequence Similarity analysis are:\nSliding Window Size,\nWindow Step Size, and the Sequence Triplet Combinations.\nThey are printed below.\n\n")
                                               print(SSAnalysisParams)
                                             }
                                             if( Step == "BlockDetection"){
                                               cat("Parameters for the detection of putative recombination blocks are:\n\nA vector containing Manual Sequence Similarity Thresholds (%),\n\nWhether or not you want HybRIDS to autodetect the sequence similarity thresholds,
                                                   \nWhether you want HybRIDS to fall back and rely on the manual thresholds should the threshold autodetection fail.
                                                   \nand finally a value by which the Standard Deviation of all sequence similarity is divded by during threshold detection (lower values = more conservative detection).\n\nThey are printed below.\n\n")
                                               print(BlockDetectionParams)
                                             }
                                             if( Step == "BlockDating"){
                                               cat("Parameters for the signifcance testing and putative recombination blocks are:\n\nA mutation rate assumed to be the average mutation rate for the sequences in the triplet,
                                                   \nA P-Value threshold - only putative recombination blocks that pass significance testing with a P-Value below the threshold are dated and output.\n\n")
                                               print(BlockDatingParams)
                                             }
                                             if( Step == "TripletGeneration"){
                                               cat("Parameters for the pre-SSAnalysis triplet generation are:\n\nThe method of triplet generation and the Sort Threshold for Method 2.
                                                   \nTo learn more about the three methods consult the documentation or HybRIDS website.\n\n")
                                               print(TripletParams)
                                             }
                                             if( Step == "Plotting"){
                                               print(PlottingParams)
                                             }
                                           },
                                         
                                         # Method for setting any parameter for any stage.
                                         setParameters =
                                           function( Step, ... ) {
                                             if( Step != "TripletGeneration" && Step != "SSAnalysis" && Step != "BlockDetection" && Step != "BlockDating" && Step != "Plotting" ){
                                               stop("You need to specify a valid analysis 'Step' to alter the paramerters of.\nThe steps are TripletGeneration, SSAnalysis, BlockDetection, BlockDating, and Plotting.")
                                             }
                                             Parameters <- list( ... )
                                             if( Step == "TripletGeneration" ){
                                               for( n in 1:length( Parameters ) ){
                                                 whichparam <- which( names( TripletParams ) == names( Parameters )[[n]])
                                                 if( class( TripletParams[[whichparam]] ) == class( Parameters[[n]] ) && length(TripletParams[[whichparam]]) == length(Parameters[[n]] ) ) {
                                                   TripletParams[[whichparam]] <<- Parameters[[n]]
                                                 } else {
                                                   warning( paste("Tried to re-assign Triplet Generation parameter ", names(TripletParams)[[whichparam]],
                                                                  " but the class of the replacement parameter or the length of the replacement parameter did not match,\nthis parameter was not changed.", sep=""))
                                                 }
                                               }
                                             }
                                             if( Step == "SSAnalysis" ) {
                                               for( n in 1:length( Parameters ) ){
                                                 whichparam <- which( names( SSAnalysisParams ) == names( Parameters )[[n]])
                                                 if( class( SSAnalysisParams[[whichparam]] ) == class( Parameters[[n]] ) && length(SSAnalysisParams[[whichparam]]) == length(Parameters[[n]] ) ) {
                                                   SSAnalysisParams[[whichparam]] <<- Parameters[[n]]
                                                 } else {
                                                   warning( paste("Tried to re-assign SSAnalysis parameter ", names(SSAnalysisParams)[[whichparam]],
                                                                  " but the class of the replacement parameter or the length of the replacement parameter did not match,\nthis parameter was not changed.", sep=""))
                                                 }
                                               } 
                                             }
                                             if( Step == "BlockDetection" ) {
                                               for( n in 1:length( Parameters ) ){
                                                 whichparam <- which( names( BlockDetectionParams ) == names( Parameters )[[n]])
                                                 if( class( BlockDetectionParams[[whichparam]] ) == class( Parameters[[n]] ) && length(BlockDetectionParams[[whichparam]]) == length(Parameters[[n]] ) ) {
                                                   BlockDetectionParams[[whichparam]] <<- Parameters[[n]]
                                                 } else {
                                                   warning( paste("Tried to re-assign Block Detection parameter ", names(BlockDetectionParams)[[whichparam]],
                                                                  " but the class of the replacement parameter or the length of the replacement parameter did not match,\nthis parameter was not changed.", sep=""))
                                                 }
                                               }
                                             }
                                             if( Step == "BlockDating" ) {
                                               for( n in 1:length( Parameters ) ){
                                                 whichparam <- which( names( BlockDatingParams ) == names( Parameters )[[n]])
                                                 if( class( BlockDatingParams[[whichparam]] ) == class( Parameters[[n]] ) && length(BlockDatingParams[[whichparam]]) == length(Parameters[[n]] ) ) {
                                                   BlockDetectionParams[[whichparam]] <<- Parameters[[n]]
                                                 } else {
                                                   warning( paste("Tried to re-assign Block Detection parameter ", names(BlockDatingParams)[[whichparam]],
                                                                  " but the class of the replacement parameter or the length of the replacement parameter did not match,\nthis parameter was not changed.", sep=""))
                                                 }
                                               }
                                             }
                                             if( Step == "Plotting" ) {
                                               for( n in 1:length( Parameters ) ){
                                                 whichparam <- which( names( PlottingParams ) == names( Parameters )[[n]])
                                                 if( names( Parameters )[[n]] == "What" && class( PlottingParams[[whichparam]] ) == class( Parameters[[n]] ) ){
                                                   PlottingParams[[whichparam]] <<- Parameters[[n]]
                                                 } else {
                                                   if( class( PlottingParams[[whichparam]] ) == class( Parameters[[n]] ) && length(PlottingParams[[whichparam]]) == length(Parameters[[n]] ) ) {
                                                     PlottingParams[[whichparam]] <<- Parameters[[n]]
                                                   } else {
                                                     warning( paste("Tried to re-assign Plotting parameter ", names(PlottingParams)[[whichparam]],
                                                                    " but the class of the replacement parameter or the length of the replacement parameter did not match,\nthis parameter was not changed.", sep=""))
                                                   }
                                                 }
                                               }
                                             }
                                           },
                                         
                                         # Method for analyzing the sequence similarity of triplets of sequences.
                                         analyzeSS = 
                                           function( Selections = "all" ) {
                                             if( !is.character( Selections ) ) stop( "option 'Selections' must be 'all' or a vector of the sequence triplets you want to use e.g. 'Seq1:Seq2:Seq3'" )
                                             if( length(Selections) == 1 && Selections == "all" ) {
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
                                         
                                         # Method to execute the putative block finds.
                                         findBlocks =
                                           function( Selections = "all" ){
                                             if( !is.character( Selections ) ) stop( "option 'Selections' must be 'all' or a vector of the sequence triplets you want to use e.g. 'Seq1:Seq2:Seq3'" )
                                             if( length(Selections) == 1 && Selections == "all" ) {
                                               if( length( Triplets ) < 2 ) {
                                                 cat( "Only one triplet to find the potential blocks in...\n" )
                                                 Triplets[[1]]$putativeBlockFind(BlockDetectionParams)
                                               } else {
                                                 cat( "Finding potential blocks in all the triplets...\n" )
                                                 for( i in 1:length( Triplets ) ) {
                                                   cat( "Now finding potential blocks for triplet", unlist(SSAnalysisParams$TripletCombinations[i], "\n") )
                                                   Triplets[[i]]$putativeBlockFind(BlockDetectionParams)
                                                 }
                                               }
                                             } else {
                                               indexTriplets( Selections )
                                               for( i in LastTripletSelection ){
                                                 cat( "Now finding potential blocks in triplet", unlist(SSAnalysisParams$TripletCombinations[i]), "\n" )
                                                 Triplets[[i]]$putativeBlockFind(BlockDetectionParams)
                                               }
                                             }
                                           },
                                         
                                         # Method to Date the blocks found.
                                         dateBlocks =
                                           function( Selections = "all" ){
                                             if( !is.character( Selections ) ) stop( "option 'Selections' must be 'all' or a vector of the sequence triplets you want to use e.g. 'Seq1:Seq2:Seq3'" )
                                             if( length(Selections) == 1 && Selections == "all" ) {
                                               if( length( Triplets ) < 2 ) {
                                                 cat( "Only one triplet to date blocks in...\n" )
                                                 Triplets[[1]]$blockDate(DNA, BlockDatingParams)
                                               } else {
                                                 cat( "Assessing and dating blocks in all the triplets...\n" )
                                                
                                                 for( i in 1:length( Triplets ) ) {
                                                   cat( "Now assessing and dating blocks for triplet", unlist(SSAnalysisParams$TripletCombinations[i], "\n"))
                                                   Triplets[[i]]$blockDate(DNA, BlockDatingParams)
                                                 }
                                               }
                                             } else {
                                               indexTriplets( Selections )
                                               for( i in LastTripletSelection ){
                                                 cat( "Now assessing and dating blocks in triplet", unlist(SSAnalysisParams$TripletCombinations[i]), "\n" )
                                                 Triplets[[i]]$blockDate(DNA, BlockDatingParams)
                                               }
                                             }
                                           },
                                         
                                         
                                         # GGplot method for HybRIDS object - activates sub-methods of triplets.
                                         plotSS =
                                           function( Selections, Combine = TRUE, ReplaceParams = TRUE, ... ) {
                                             if( !is.character( Selections ) ) stop( "option 'Selections' must be a vector of the sequence triplets you want to use e.g. 'Seq1:Seq2:Seq3'" )
                                             oldParameters <- PlottingParams
                                             newParameters <- list( ... )
                                             if(length(newParameters) > 0){
                                               for( n in 1:length(newParameters) ){
                                                 whichparam <- which(names(PlottingParams) == names(newParameters)[[n]])
                                                 if(class(PlottingParams[[whichparam]]) == class(newParameters[[n]])){
                                                   PlottingParams[[whichparam]] <<- newParameters[[n]]
                                                 } else {
                                                   warning( paste("Tried to re-assign plotting parameter ", names(PlottingParams)[[whichparam]],
                                                                  " but the class of the replacement parameter did not match, this parameter was not changed.", sep=""))
                                                 }
                                               }
                                             }
                                             for( i in Selections ) {
                                               cat("Selection", i)
                                               if( length( unlist( strsplit(i, ":") ) ) == 3 ) {
                                                 indexTriplets( i )
                                                 if( "Lines" %in% PlottingParams$What && !"Bars" %in% PlottingParams$What ) {
                                                   outplot <- Triplets[[LastTripletSelection]]$plotLines( PlottingParams )
                                                 }
                                                 if( !"Lines" %in% PlottingParams$What && "Bars" %in% PlottingParams$What ) {
                                                   outplot <- Triplets[[LastTripletSelection]]$plotBars( parameters = PlottingParams )
                                                 }
                                                 if( "Lines" %in% PlottingParams$What && "Bars" %in% PlottingParams$What ) {
                                                   outplot <- arrangeGrob( Triplets[[LastTripletSelection]]$plotBars( parameters = PlottingParams ),
                                                                Triplets[[LastTripletSelection]]$plotLines( PlottingParams ),
                                                                ncol = 1 )
                                                 }
                                                 if(ReplaceParams == FALSE){
                                                   PlottingParams <<- oldParameters
                                                 }
                                                 return(outplot)
                                               } else {
                                                 if( length( unlist( strsplit( i, ":" ) ) ) == 2 && Combine == TRUE ) {
                                                   indexTriplets( i )
                                                   if( "Lines" %in% PlottingParams$What ){
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
                                                       geom_line( aes( colour = TripletSet ), show_guide = PlottingParams$Legends, size = 0.8 ) +
                                                       ylab( "% Sequence Similarity" ) +
                                                       xlab( "Base Position" )
                                                     outplotLines <- applyPlottingParams( outplotLines, PlottingParams, title = paste("Sequence similarity for sequence pair ", i, " in all triplets in which it occurs", sep="" ) )
                                                   }
                                                   if( "Bars" %in% PlottingParams$What ){
                                                     bars <- lapply( Triplets[LastTripletSelection], function(x) x$plotBars( exportDat = T, PlottingParams ) )
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
                                                       ylab( "Triplet Set" ) +
                                                       scale_x_continuous( breaks = c(seq( from = 1, to = PlottingParams$MosaicScale, by = PlottingParams$MosaicScale / 10 ), PlottingParams$MosaicScale), labels = c(bpX[seq( from = 1, to = PlottingParams$MosaicScale, by = PlottingParams$MosaicScale / 10 )], max(bpX)) ) +
                                                       scale_y_discrete( labels = as.character(yaxislab) ) +
                                                       scale_fill_gradient2(high="red", low="blue", midpoint=33.3)
                                                     outplotBars <- applyPlottingParams(outplotBars, PlottingParams, title = paste("Sequence similarity for sequence pair ", i, " in all triplets in which it occurs", sep="" ) )
                                                   }
                                                   if( "Lines" %in% PlottingParams$What && !"Bars" %in% PlottingParams$What ) {
                                                     if(ReplaceParams == FALSE){
                                                       PlottingParams <<- oldParameters
                                                     }
                                                     return(outplotLines)
                                                   } else {
                                                     if( !"Lines" %in% PlottingParams$What && "Bars" %in% PlottingParams$What ) {
                                                       if(ReplaceParams == FALSE){
                                                         PlottingParams <<- oldParameters
                                                       }
                                                       return(outplotBars)
                                                     } else {
                                                       if( "Lines" %in% PlottingParams$What && "Bars" %in% PlottingParams$What && Combine == TRUE ) {
                                                         if( PlottingParams$CombinedTitle == TRUE ){
                                                           if(ReplaceParams == FALSE){
                                                             PlottingParams <<- oldParameters
                                                           }
                                                           return( arrangeGrob( textGrob( paste( "Sequence similarity for sequence pair ", i, " in all triplets in which it occurs", sep="" ), x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                                                                                         just = "centre" ), outplotBars, outplotLines, ncol = 1) )
                                                         } else {
                                                           output <- arrangeGrob( outplotBars, outplotLines, ncol = 1 )
                                                           if(ReplaceParams == FALSE){
                                                             PlottingParams <<- oldParameters
                                                           }
                                                           return( output )
                                                         }
                                                       } else {
                                                         if(ReplaceParams == FALSE){
                                                           PlottingParams <<- oldParameters
                                                         }
                                                         return( list( Barplot = outplotBars, Linesplot = outplotLines) )
                                                       }
                                                     }
                                                   }
                                                 } else {
                                                   if(length( unlist( strsplit( i, ":" ) ) ) == 2 && Combine == FALSE){
                                                     
                                                   }
                                                 }
                                               }
                                             }
                                             },
                                         
                                         # Method to put the data from detected blocks in triplets into a data format.
                                         tabulateDetectedBlocks =
                                           function( Selection, OneTable = FALSE ) {
                                             outputTables <- list()
                                             len <- length( unlist( lapply( Selection, function(x) indexTriplets( x, output = TRUE ) ) ) )
                                             ind <- unlist( lapply( Selection, function(x) indexTriplets( x, output = TRUE ) ) )
                                             tables <- lapply( Triplets[ind], function(x) x$tabulateBlocks() )
                                                                                
                                             
                                             return(tables)
                                           },
                                                        
                                         # Method for indexing triplets.
                                         indexTriplets =
                                           function( selections, output = F ) {
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
                                             if(output == T){
                                               return(LastTripletSelection)
                                             }
                                           },
                                         
                                         # Show method.
                                         show = function() {
                                           cat("HybRIDS object - Analysis of ",length(DNA$SequenceNames)," aligned sequences.\n\nDNA Alignment:\n--------------\nFull Sequence File Location: ", DNA$FullSequenceFile,
                                               "\nInformative Sequence File Location: ", DNA$InformativeSequenceFile) 
                                           if(length(DNA$SequenceNames) > 0){
                                             cat("\nFull Length: ",DNA$SequenceLength,
                                                 "\nInformative Length: ",DNA$InformativeLength,"\nSequence names: ",DNA$SequenceNames,"\n\n")
                                           } else {
                                             cat("\n\n A DNA sequence alignment file has not yet been loaded into the HybRIDS object.\n\n")
                                           }
                                           cat("Triplet Generation Parameters:\n------------------------------\nTriplet Generation Method: ",TripletParams$Method,
                                               "\nThreshold for method number 2: ",TripletParams$SortThreshold,
                                               "\n\nSequence Similarity Analysis Parameters:\n----------------------------------------\n",
                                               "Sliding Window Size: ",SSAnalysisParams$WindowSize,"\nSliding Window Step Size: ",SSAnalysisParams$StepSize,
                                               "\n\nBlock Detection Parameters: \n---------------------------\nManual Thresholds: ",BlockDetectionParams$ManualThresholds,sep="")
                                           if(BlockDetectionParams$AutoThresholds == TRUE){
                                             cat("\nHybRIDS will attempt automatic detection of SS Thresholds for putative block searches.")
                                             if(BlockDetectionParams$ManualFallback == TRUE){
                                               cat("\nHybRIDS will fall back on user specified manual thresholds, should the autodetection fail.")
                                             } else {
                                               cat("\nHybRIDS will not fall back on user specified manual thresholds, should the autodetection fail.")
                                             }
                                           } else {
                                             cat("\nHybRIDS will not attempt automatic detection of SS Thresholds for putative block searches.\nOnly the manually specified thresholds will be used.")
                                           }
                                           cat("\n\nBlock Dating Parameters:\n------------------------\nAssumed mutation rate: ",BlockDatingParams$MutationRate,
                                               "\nP-Value for acceptance of putative blocks: ",BlockDatingParams$PValue,sep="")
                                           if( length(Triplets) < 1 ){
                                             cat("\n\nNo Triplets have been generated with the method makeTripletCombos yet.")
                                           } else {
                                             cat("\n\n",length(Triplets),"Triplet(s) have been generated for analysis.")
                                           }  
                                         }
                                         



                                          
                          )
                        
                        )














                          