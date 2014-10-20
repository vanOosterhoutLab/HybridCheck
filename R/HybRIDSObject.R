#' @docType package
#' ...
#' @useDynLib HybRIDS
#' @import ape ggplot2 grid gridExtra png
NULL

#' A Reference Class for manageing a HybRIDS analysis.
#' @name HybRIDS
#' @import methods
#' @export HybRIDS
#' @exportClass HybRIDS
#' @field DNA A HybRIDSdna reference object.
HybRIDS <- setRefClass("HybRIDS",
                        
                        fields = list( 
                          DNA = "ANY",
                          ssAnalysisSettings = "ANY",
                          BlockDetectionParams = "list",
                          BlockDatingParams = "list",
                          LastTripletSelection = "numeric",
                          PlottingParams = "list",
                          InGUI = "logical",
                          triplets = "ANY",
                          userBlocks = "ANY",
                          filesDirectory = "character"
                          ),
                        
                         methods = list(initialize = 
                                           function(dnafile=NULL, inGUI=FALSE){
                                             "Create HybRIDS object with default values for fields. The path to the FASTA file can be provided."
                                             filesDirectory <<- tempdir()
                                             InGUI <<- inGUI
                                             
                                             # Initiate settings for SSAnalysis scans.
                                             ssAnalysisSettings <<- SSAnalysisSettings$new()
                                             
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
                                             
                                             # Initialize the user blocks object.
                                             userBlocks <<- UserBlocks$new()
                                             
                                             # Initialize the HybRIDSseq object.
                                             DNA <<- HybRIDSseq$new()
                                             
                                             triplets <<- Triplets$new() 
                                             
                                             # Make sure the input DNA file is loaded into the HybRIDSseq object.
                                             if(!is.null(dnafile)){
                                               DNA$InputDNA(dnafile)
                                               userBlocks$initializePairsFromDNA(DNA)
                                             } else {
                                               stop("You haven't provided a path to a DNA file or the name of an object of class DNAbin")
                                             }
                                           },
                                         
                                         # Method for generating the Triplet Combinations...
                                         generateComparrisons =
                                           function(...){
                                             triplets$decideTriplets(DNA, list(...))
                                           },
                                         
                                         # Method for setting any parameter for any stage.
                                         setParameters =
                                           function(Step = NULL, ...){
                                             if(!any(Step == c("TripletGeneration", "SSAnalysis", "BlockDetection", "BlockDating", "Plotting"))){
                                               stop("You need to specify a valid analysis 'Step' to alter the paramerters of.\nThe steps are TripletGeneration, SSAnalysis, BlockDetection, BlockDating, and Plotting.")
                                             }
                                             Parameters <- list(...)
                                             if(Step == "TripletGeneration"){
                                               triplets$comparrisonSettings$setSettings(...)
                                             }
                                             if(Step == "SSAnalysis") {
                                               ssAnalysisSettings$setSettings(...) 
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
                                             DNA$enforceDNA()
                                             if(!is.character(Selections)) stop("option 'Selections' must be 'ALL' or a vector of the sequence triplets you want to use e.g. 'Seq1:Seq2:Seq3'")
                                             if(!comparrisonSettings$hasMultipleCombinations()){
                                                message("Only one triplet to analyze the sequence similarity of...")
                                                seq.similarity(DNA, Triplets[[1]], ssAnalysisSettings)
                                              } else {
                                                indexTriplets(Selections)
                                                for(i in LastTripletSelection){
                                                  message("Now analysing sequence similarity of triplet ", names(Triplets)[i])
                                                  suppressMessages(seq.similarity(DNA, Triplets[[i]], ssAnalysisSettings))
                                               } 
                                             }
                                             message("Finished Sequence Similarity Analysis.")
                                           },
                                         
                                         # Method to execute the putative block finds.
                                         findBlocks =
                                           function(Selections = "ALL"){
                                             if(!is.character(Selections)) stop( "option 'Selections' must be 'ALL' or a vector of the sequence triplets you want to use e.g. 'Seq1:Seq2:Seq3'" )
                                             if(length(Triplets) < 2){
                                               message("Only one triplet to find the potential blocks in...")
                                               Triplets[[1]]$putativeBlockFind(BlockDetectionParams)
                                             } else {
                                               indexTriplets(Selections)
                                               for(i in LastTripletSelection){
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
                                               indexTriplets(Selections)
                                               for(i in LastTripletSelection){
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
                                           function(selections, output=F) {
                                             if(!is.character(selections)) stop("option 'which' must be a vector of the sequence triplets you want to use e.g. 'Seq1:Seq2:Seq3'")
                                             if(any(selections == "ALL")){
                                               message("One of the selections is 'ALL', any other selections provided are redundant...")
                                               LastTripletSelection <<- 1:length(Triplets)
                                               if(output==T){
                                                 return(LastTripletSelection)
                                               }
                                             } else {
                                               processedSelections <- strsplit(selections, split=":")
                                               threes <- processedSelections[which(lapply(processedSelections, function(x) length(x)) == 3)]
                                               twos <- processedSelections[which(lapply(processedSelections, function(x) length(x)) == 2)]
                                               ones <- processedSelections[which(lapply(processedSelections, function(x) length(x)) == 1)]
                                               threes <- lapply(threes, function(x) which(DNA$SequenceNames %in% x))
                                               twos <- lapply(twos, function(x) which(DNA$SequenceNames %in% x))
                                               ones <- lapply(ones, function(x) which(DNA$SequenceNames %in% x))
                                               threes <- which(SSAnalysisParams$TripletCombinations %in% threes)
                                               if(length(twos) > 0){
                                                 twos <- which(unlist(lapply(twos, function(y) lapply(SSAnalysisParams$TripletCombinations, function(x) all(y %in% x)))))
                                               } else {
                                                 twos <- c()
                                               }
                                               if(length(ones) > 0){
                                                 ones <- which(unlist(lapply(ones, function(y) lapply(SSAnalysisParams$TripletCombinations, function(x) any(y %in% x)))))
                                               } else {
                                                 ones <- c()
                                               }
                                               # Now let's get rid of redunancies and assign the selection to LastTripletSelection.
                                               LastTripletSelection <<- unique(c(threes, twos, ones))
                                               if(output == T){
                                                 return(LastTripletSelection)
                                               } 
                                             }
                                           },
                                         
                                         # Show method.
                                         show = function() {
                                           cat("HybRIDS object:\n\n")
                                           DNA$show()
                                           cat("\n\n")
                                           triplets$show()
                                           cat("\n\n")
                                            
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
                                         
                                         addUserBlock = function(selection, firstbp, lastbp){
                                           userBlocks$addBlock(firstbp, lastbp, selection)
                                         },
                                        
                                         clearUserBlocks = function(selection){
                                           userBlocks$blankBlocks(selection)
                                         },
                                        
                                         dateUserBlocks = function(){
                                           userBlocks$dateBlocks(DNA, BlockDatingParameters)
                                         },
                                         
                                         recombinationExtent = function(Selection = "ALL"){
                                           return((sum(tabulateDetectedBlocks(Selection)$Approximate_Length_BP)/DNA$SequenceLength)*100)
                                         }
                          )
                        )

