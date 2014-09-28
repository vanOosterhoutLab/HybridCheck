#' @docType package
#' ...
#' @import ape ggplot2 grid gridExtra png
NULL

#' @title HybRIDS reference class
#' @name HybRIDS
#' @description 
#' The HybRIDS reference class is the main class that is used
#' @export
HybRIDS <- setRefClass("HybRIDS",
                        
                        fields = list( 
                          DNA = "ANY",
                          TripletParams = "list",
                          SSAnalysisParams = "list",
                          BlockDetectionParams = "list",
                          BlockDatingParams = "list",
                          LastTripletSelection = "numeric",
                          PlottingParams = "list",
                          InGUI = "logical",
                          Triplets = "list",
                          UserBlocks = "list"
                          ),
                        
                         methods = list(initialize = 
                                           function(dnafile=NULL, formatForce=NULL, storageOpt="default", inGUI=FALSE){
                                             # Initiate settings for triplet generation.
                                             TripletParams <<- list(
                                               Method = 1,
                                               DistanceThreshold = 0.01,
                                               PartitionStrictness = 2,
                                               Refine = FALSE)
                                             # Initiate settings for sliding window scans.
                                             SSAnalysisParams <<- list(
                                               WindowSize = 100,
                                               StepSize = 1,
                                               TripletCombinations = list())
                                             #length(BlockDetectionParams) <<- 4
                                             # Initiate block detection parameters.
                                             BlockDetectionParams <<- list(ManualThresholds = c(90), AutoThresholds = TRUE, ManualFallback = TRUE, SDstringency = 2)
                                             # Initiate block dating parameters.
                                             BlockDatingParams <<- list(MutationRate = 10e-08, PValue = 0.005, BonfCorrection = TRUE, DateAnyway = FALSE)
                                             # Initiate settings for plotting triplet
                                             PlottingParams <<- list(What = c("Bars", "Lines"), PlotTitle = TRUE, CombinedTitle = FALSE, 
                                                                      TitleSize = 14, TitleFace="bold", TitleColour = "black", XLabels = TRUE, YLabels = TRUE,
                                                                      XTitle = TRUE, XTitleFontSize = 12, XTitleColour = "black", XLabelSize = 10, XLabelColour = "black",
                                                                      YTitle = TRUE, YTitleFontSize = 12, YTitleColour = "black", YLabelSize = 10, YLabelColour = "black",
                                                                      Legends = TRUE, LegendFontSize = 12, MosaicScale = 500)
                                             InGUI <<- inGUI
                                             DNA <<- HybRIDSseq$new()
                                             if(!is.null(dnafile)){
                                               DNA$InputDNA(dnafile, formatForce)
                                               pairs <- unlist(lapply(combn(unique(DNA$SequenceNames),2, simplify=F), function(x) paste(x[1], x[2], sep=":")))
                                               lapply(pairs, function(x) UserBlocks[[x]] <<- data.frame(FirstBP=as.numeric(), LastBP=as.numeric(), ApproxBpLength=as.numeric()))
                                             } else {
                                               stop("You haven't provided a path to a DNA file or the name of an object of class DNAbin")
                                             }
                                           },
                                         
                                         # Method for generating the Triplet Combinations...
                                         makeTripletCombos =
                                           function(...){
                                             # Input checking - check DNA sequences have been read in, and 
                                             if(!DNA$hasDNA()) stop("No DNA data is loaded into this HybRIDS object")
                                             if(nrow(DNA$InformativeSequence) < 3){
                                               if(InGUI == TRUE) gmessage("Not enough sequences to make any triplets", icon="error")
                                               stop("Not enough sequences to make any triplets, most likely the removal of duplicate sequences has resulted in too few sequences.")
                                             }
                                             if(!TripletParams$Refine){
                                               SSAnalysisParams$TripletCombinations <<- combn(c(1:nrow(DNA$InformativeSequence)), 3, simplify=FALSE) 
                                             }
                                             rejects <- c()
                                             # Triplet Generation Method 1 is simply to include all possible triplets.
                                             # In which case the below code is skipped...
                                             if(TripletParams$Method > 1 && length(SSAnalysisParams$TripletCombinations) > 1){
                                               ranges <- list(...)
                                               if(TripletParams$Method == 2){
                                                 # Use the partition method to generate triplets to check for recombination between partitions.
                                                 if(length(ranges) <= 1){
                                                   stop("Error: You need to provide more than one valid partition.")
                                                 }
                                                 if(any(unlist(ranges) > nrow(DNA$InformativeSequence))){
                                                   stop("Error: Provided sequence numbers in the partitions, higher than the actual number of sequences in alignment.")
                                                 }
                                                 message("Generating triplets to find recombination in sequences, between partitions.")
                                                 rejects <- unlist(lapply(ranges, function(x){
                                                   which(unlist(lapply(SSAnalysisParams$TripletCombinations, function(y){
                                                     length(which(x %in% y)) > TripletParams$PartitionStrictness
                                                   })))
                                                 }))
                                               }
                                               # Use the method whereby distance information is used to reject pairs which would likeley be pointless.
                                               if(TripletParams$Method == 3 || TripletParams$Method == 4){                                                    
                                                 distances <- dist.dna(as.DNAbin(DNA$FullSequence), model = "raw")
                                                 seqpairs <- combn(c(1:nrow(DNA$InformativeSequence)), 2, simplify=FALSE)
                                                 if(TripletParams$Method == 3){
                                                   # Reject distances that are below a given threshold.
                                                   if(length(ranges) > 1 || class(ranges[[1]])){
                                                     stop("Error: One integer value should be provided as a distance threshold.")
                                                   }
                                                   rejectiondistances <- seqpairs[which(distances < TripletParams$DistanceThreshold)]
                                                 } else {
                                                   distances_density <- density(distances)
                                                   Lows <- cbind(distances_density$x[which(diff(sign(diff(distances_density$y))) == 2)], distances_density$y[which(diff(sign(diff(distances_density$y))) == 2)])
                                                   Lowest <- Lows[which(Lows[,1] == min(Lows[,1])),]
                                                   rejectiondistances <- seqpairs[which(distances < Lowest[1])]
                                                 }
                                                 for(i in 1:length(SSAnalysisParams$TripletCombinations)){
                                                   for(n in 1:length(rejectiondistances)){
                                                     if(all(rejectiondistances[[n]] %in% SSAnalysisParams$TripletCombinations[[i]])){
                                                       rejects <- append(rejects, i)
                                                       break
                                                     }
                                                   }
                                                 }
                                               }
                                             }
                                             if(!is.null(rejects) && length(rejects) > 0){
                                               SSAnalysisParams$TripletCombinations <<- SSAnalysisParams$TripletCombinations[-unlist(rejects)]
                                             }
                                             Triplets <<- lapply(SSAnalysisParams$TripletCombinations, function(x) HybRIDStriplet$new(sequencenumbers = x, sequences = c(DNA$SequenceNames[x[1]], DNA$SequenceNames[x[2]], DNA$SequenceNames[x[3]]), fullseqlength = DNA$SequenceLength))
                                             names(Triplets) <<- unlist(lapply(SSAnalysisParams$TripletCombinations, function(x) paste(DNA$SequenceNames[x[1]], DNA$SequenceNames[x[2]], DNA$SequenceNames[x[3]], sep = ":")))
                                           },
                                         
                                         # Method for displaying parameters.
                                         showParameters = 
                                           function(Step = NULL){
                                             if( Step!= "TripletGeneration" && Step != "SSAnalysis" && Step != "BlockDetection" && Step != "BlockDating" && Step != "Plotting"){
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
                                           function(Step = NULL, ...){
                                             if(Step != "TripletGeneration" && Step != "SSAnalysis" && Step != "BlockDetection" && Step != "BlockDating" && Step != "Plotting" ){
                                               stop("You need to specify a valid analysis 'Step' to alter the paramerters of.\nThe steps are TripletGeneration, SSAnalysis, BlockDetection, BlockDating, and Plotting.")
                                             }
                                             Parameters <- list(...)
                                             if(Step == "TripletGeneration"){
                                               for(n in 1:length(Parameters)){
                                                 whichparam <- which(names(TripletParams) == names(Parameters)[[n]])
                                                 if(class(TripletParams[[whichparam]]) == class(Parameters[[n]]) && length(TripletParams[[whichparam]]) == length(Parameters[[n]])){
                                                   TripletParams[[whichparam]] <<- Parameters[[n]]
                                                 } else {
                                                   warning(paste("Tried to re-assign Triplet Generation parameter ", names(TripletParams)[[whichparam]],
                                                                  " but the class of the replacement parameter or the length of the replacement parameter did not match,\nthis parameter was not changed.", sep=""))
                                                 }
                                               }
                                             }
                                             if(Step == "SSAnalysis") {
                                               for(n in 1:length(Parameters)){
                                                 whichparam <- which(names(SSAnalysisParams) == names(Parameters)[[n]])
                                                 if(class(SSAnalysisParams[[whichparam]]) == class(Parameters[[n]]) && length(SSAnalysisParams[[whichparam]]) == length(Parameters[[n]])){
                                                   SSAnalysisParams[[whichparam]] <<- Parameters[[n]]
                                                 } else {
                                                   warning( paste("Tried to re-assign SSAnalysis parameter ", names(SSAnalysisParams)[[whichparam]],
                                                                  " but the class of the replacement parameter or the length of the replacement parameter did not match,\nthis parameter was not changed.", sep=""))
                                                 }
                                               } 
                                             }
                                             if(Step == "BlockDetection") {
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
                                             if(Step == "BlockDating") {
                                               for(n in 1:length(Parameters)){
                                                 whichparam <- which(names(BlockDatingParams) == names(Parameters)[[n]])
                                                 if(class(BlockDatingParams[[whichparam]]) == class(Parameters[[n]]) && length(BlockDatingParams[[whichparam]]) == length(Parameters[[n]])){
                                                   BlockDatingParams[[whichparam]] <<- Parameters[[n]]
                                                 } else {
                                                   warning(paste("Tried to re-assign Block Detection parameter ", names(BlockDatingParams)[[whichparam]],
                                                                  " but the class of the replacement parameter or the length of the replacement parameter did not match,\nthis parameter was not changed.", sep=""))
                                                 }
                                               }
                                             }
                                             if(Step == "Plotting") {
                                               for(n in 1:length(Parameters)){
                                                 whichparam <- which(names(PlottingParams) == names(Parameters)[[n]])
                                                 if(names(Parameters)[[n]] == "What" && class(PlottingParams[[whichparam]]) == class(Parameters[[n]])){
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
                                           function(Selections = "ALL"){
                                             if(!is.character(Selections)) stop("option 'Selections' must be 'ALL' or a vector of the sequence triplets you want to use e.g. 'Seq1:Seq2:Seq3'")
                                             if(length(SSAnalysisParams$TripletCombinations) < 2){
                                                message("Only one triplet to analyze the sequence similarity of...")
                                                seq.similarity(DNA$InformativeSequence, Triplets[[1]], SSAnalysisParams$WindowSize, SSAnalysisParams$StepSize, DNA$SequenceLength, DNA$InformativeBp)
                                              } else {
                                                indexTriplets(Selections)
                                                for(i in LastTripletSelection){
                                                  message("Now analysing sequence similarity of triplet ", paste(unlist(SSAnalysisParams$TripletCombinations[i]), collapse=":"))
                                                  suppressMessages(seq.similarity(DNA$InformativeSequence[unlist(SSAnalysisParams$TripletCombinations[[i]]),], Triplets[[i]], SSAnalysisParams$WindowSize, SSAnalysisParams$StepSize, DNA$SequenceLength, DNA$InformativeBp))
                                               } 
                                             }
                                             message("Finished Sequence Similarity Analysis.")
                                           },
                                         
                                         # Method to execute the putative block finds.
                                         findBlocks =
                                           function( Selections = "ALL" ){
                                             if( !is.character( Selections ) ) stop( "option 'Selections' must be 'ALL' or a vector of the sequence triplets you want to use e.g. 'Seq1:Seq2:Seq3'" )
                                             if( length( Triplets ) < 2 ) {
                                               message("Only one triplet to find the potential blocks in...")
                                               Triplets[[1]]$putativeBlockFind(BlockDetectionParams)
                                             } else {
                                               indexTriplets( Selections )
                                               for( i in LastTripletSelection ){
                                                 message("Now finding potential blocks for triplet ", paste(unlist(SSAnalysisParams$TripletCombinations[i]), collapse=":"))
                                                 suppressMessages(Triplets[[i]]$putativeBlockFind(BlockDetectionParams))
                                               }
                                             }
                                           },
                                         
                                         # Method to Date the blocks found.
                                         dateBlocks =
                                           function(Selections = "ALL"){
                                             if(!is.character(Selections)) stop("option 'Selections' must be 'ALL' or a vector of the sequence triplets you want to use e.g. 'Seq1:Seq2:Seq3'")
                                             if(length(Triplets) < 2){
                                               message("Only one triplet to date blocks in...")
                                               Triplets[[1]]$blockDate(DNA, BlockDatingParams)
                                             } else {
                                               indexTriplets( Selections )
                                               for( i in LastTripletSelection ){
                                                 message("Now assessing and dating blocks for triplet ", paste(unlist(SSAnalysisParams$TripletCombinations[i]), collapse=":"))
                                                 suppressMessages(Triplets[[i]]$blockDate(DNA, BlockDatingParams))
                                               }
                                             }
                                           },
                                         
                                         
                                         # GGplot method for HybRIDS object - activates sub-methods of triplets.
                                         plotTriplets =
                                           function(Selections = "ALL", Combine = TRUE, ReplaceParams = TRUE, ...){
                                             if(!is.character(Selections)) stop("option 'Selections' must be a vector of the sequence triplets you want to use e.g. 'Seq1:Seq2:Seq3'")
                                             oldParameters <- PlottingParams
                                             newParameters <- list(...)
                                             if(length(newParameters) > 0){
                                               for(n in 1:length(newParameters)){
                                                 whichparam <- which(names(PlottingParams) == names(newParameters)[[n]])
                                                 if(class(PlottingParams[[whichparam]]) == class(newParameters[[n]])){
                                                   PlottingParams[[whichparam]] <<- newParameters[[n]]
                                                 } else {
                                                   warning(paste("Tried to re-assign plotting parameter ", names(PlottingParams)[[whichparam]],
                                                                  " but the class of the replacement parameter did not match, this parameter was not changed.", sep=""))
                                                 }
                                               }
                                             }
                                             indexTriplets(Selections)
                                             if("Lines" %in% PlottingParams$What && "Bars" %in% PlottingParams$What && Combine == TRUE){
                                               outplot <- lapply(LastTripletSelection, function(i){Triplets[[i]]$combineLinesAndBars(PlottingParams)})
                                             }
                                             if("Lines" %in% PlottingParams$What && "Bars" %in% PlottingParams$What && Combine == FALSE ){
                                               outplot <- lapply(LastTripletSelection, function(i){list(bars = Triplets[[i]]$plotBars(parameters = PlottingParams), lines = Triplets[[i]]$plotLines(PlottingParams))})  
                                             }
                                             if("Lines" %in% PlottingParams$What && !"Bars" %in% PlottingParams$What){
                                               outplot <- lapply(LastTripletSelection, function(i){Triplets[[i]]$plotLines(PlottingParams)})
                                             }
                                             if(!"Lines" %in% PlottingParams$What && "Bars" %in% PlottingParams$What){
                                               outplot <- lapply(LastTripletSelection, function(i){Triplets[[i]]$plotBars(parameters = PlottingParams)})
                                             }
                                             if(length(outplot) == 1){
                                               outplot <- outplot[[1]]
                                             }
                                             if(ReplaceParams == FALSE){
                                               PlottingParams <<- oldParameters
                                             }
                                             return(outplot)
                                             },
                                         
                                         # Method to put the data from detected blocks in triplets into a data format.
                                         tabulateDetectedBlocks =
                                           function(Selection = "ALL", OneTable = FALSE, Neat = TRUE) {
                                             outputTables <- list()
                                             ind <- unlist(lapply(Selection, function(x) indexTriplets(x, output = TRUE)))
                                             tables <- lapply(Triplets[ind], function(x) x$tabulateBlocks())
                                             tripletlabels <- unlist(lapply(1:length(tables), function(i) rep(names(tables)[[i]], nrow(tables[[i]]))))                                
                                             tables <- do.call(rbind, tables)
                                             tables["Triplet"] <- tripletlabels
                                             output <- data.frame(tables$SequencePair, tables$SequenceSimilarityThreshold, tables$Triplet, tables$Length,
                                                                  tables$First, tables$Last, tables$FirstBP, tables$LastBP, tables$ApproxBpLength, tables$SNPnum, tables$fiveAge, tables$fiftyAge,
                                                                  tables$ninetyfiveAge, tables$PValue, tables$PThresh, tables$MeanAge, tables$CorrectedSNPs) 
                                             if(Neat == TRUE) {
                                               output <- output[,-c(4,5,6)]
                                               names(output) <- c("Sequence_Pair","Sequence_Similarity_Threshold","Triplet","First_BP_Position","Last_BP_Position","Approximate_Length_BP","Number_of_SNPs","p=0.05_Age","p=0.5_Age","p=0.95_Age","P_Value", "P_Thresh", "Mean_Age", "Corrected_Number_of_SNPs")
                                             }
                                             class(output) <- c(class(output), "HybRIDStable")   
                                             return(output)
                                           },
                                                        
                                         # Method for indexing triplets.
                                         indexTriplets =
                                           function( selections, output = F ) {
                                             if( !is.character( selections ) ) stop( "option 'which' must be a vector of the sequence triplets you want to use e.g. 'Seq1:Seq2:Seq3'" )
                                             if(any(selections == "ALL")){
                                               message("One of the selections is 'ALL', any other selections provided are redundant...")
                                               LastTripletSelection <<- 1:length(Triplets)
                                               if(output==T){
                                                 return(LastTripletSelection)
                                               }
                                             } else {
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
                                             }
                                           },
                                         
                                         # Show method.
                                         show = function() {
                                           cat("HybRIDS object - Analysis of ",length(DNA$SequenceNames)," aligned sequences.\n\nDNA Alignment:\n--------------") 
                                           if(length(DNA$SequenceNames) > 0){
                                             cat("\nFull Length: ",DNA$SequenceLength,
                                                 "\nInformative Length: ",DNA$InformativeLength,"\nSequence names: ",DNA$SequenceNames,"\n\n")
                                           } else {
                                             cat("\n\n A DNA sequence alignment file has not yet been loaded into the HybRIDS object.\n\n")
                                           }
                                           cat("Triplet Generation Parameters:\n------------------------------\nTriplet Generation Method: ",TripletParams$Method,
                                               "\nDistance threshold for method number 3: ",TripletParams$SortThreshold,
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
                                         },
                                         
                                         userBlockAdd = function(first, last, pair) {
                                           selections <- unlist(strsplit(pair,":"))
                                           if(length(selections) != 2){ stop("You must specify two sequences, between which your recombination event occured.") }
                                           options <- strsplit(names(UserBlocks),":")
                                           index <- which(unlist(lapply(lapply(options, function(x) selections %in% x), function(y) all(y))))
                                           if(length(index) != 1){ stop("Something has gone wrong indexing pairs in triplets - this scenario should not happen, the index of more than or less than one pair should not be possible, contact package maintainer.")}
                                           bplength <- abs(last-first)+1
                                           UserBlocks[[index]] <<- rbind(UserBlocks[[index]], c(first, last, bplength))
                                           names(UserBlocks[[index]]) <<- c("FirstBP", "LastBP", "ApproxBpLength")
                                         },
                                         userBlockBlank = function(selection){
                                           selections <- unlist(strsplit(pair,":"))
                                           if(length(selections) != 2){ stop("You must specify two sequences, between which your recombination event occured.") }
                                           options <- strsplit(names(UserBlocks),":")
                                           index <- which(unlist(lapply(lapply(options, function(x) selections %in% x), function(y) all(y))))
                                           if(length(index) != 1){ stop("Something has gone wrong indexing pairs in triplets - this scenario should not happen, the index of more than or less than one pair should not be possible, contact package maintainer.")}
                                           UserBlocks[[index]] <<- data.frame(FirstBP=as.numeric(), LastBP=as.numeric(), ApproxBpLength=as.numeric())
                                         },
                                         dateCustomBlocks = function(){
                                           for(i in 1:length(UserBlocks)){
                                             if(nrow(UserBlocks[[i]]) > 0){
                                                Pair <- which(DNA$SequenceNames %in% unlist(strsplit(names(UserBlocks)[i], ":")))
                                                dated <- date.blocks(UserBlocks[[i]], DNA, BlockDatingParams$MutationRate, Pair, BlockDatingParams$PValue, BlockDatingParams$BonfCorrection, BlockDatingParams$DateAnyway)
                                                UserBlocks[[i]] <<- cbind(UserBlocks[[i]], dated)
                                           }
                                         }
                                         },
                                         
                                         recombinationExtent = function(Selection = "ALL"){
                                           return((sum(tabulateDetectedBlocks(Selection)$Approximate_Length_BP)/DNA$SequenceLength)*100)
                                         }
                          )
                        )
