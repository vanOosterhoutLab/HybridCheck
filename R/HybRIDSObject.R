#' @docType package
#' ...
#' @useDynLib HybRIDS
#' @import ape ggplot2 grid gridExtra png
NULL

#' A Reference Class for managing a HybRIDS analysis.
#' @name HybRIDS
#' @description 
#' @import methods
#' @export HybRIDS
#' @exportClass HybRIDS
#' @field DNA A HybRIDSdna reference object.
HybRIDS <- setRefClass("HybRIDS",
                        
                        fields = list( 
                          DNA = "ANY",
                          comparrisonSettings = "ANY",
                          ssAnalysisSettings = "ANY",
                          blockDetectionSettings = "ANY",
                          BlockDatingParams = "list",
                          plottingSettings = "ANY",
                          triplets = "ANY",
                          userBlocks = "ANY",
                          filesDirectory = "character"
                          ),
                        
                         methods = list(initialize = 
                                           function(dnafile=NULL){
                                             "Create HybRIDS object with default values for fields. The path to the FASTA file can be provided."
                                             filesDirectory <<- tempdir()
                                             
                                             # Initiate settings for SSAnalysis scans.
                                             ssAnalysisSettings <<- SSAnalysisSettings$new()
                                             
                                             # Initiate block detection parameters.
                                             blockDetectionSettings <<- BlockDetectionSettings$new()
                                             
                                             # Initiate block dating parameters.
                                             BlockDatingParams <<- list(MutationRate = 10e-08, PValue = 0.005, BonfCorrection = TRUE, DateAnyway = FALSE)
                                             
                                             # Initiate settings for plotting triplet
                                             plottingSettings <<- PlottingSettings$new()
                                             
                                             # Initialize the user blocks object.
                                             userBlocks <<- UserBlocks$new()
                                             
                                             # Initialize the HybRIDSseq object.
                                             DNA <<- HybRIDSseq$new()
                                             
                                             triplets <<- Triplets$new() 
                                             
                                             # Make sure the input DNA file is loaded into the HybRIDSseq object.
                                             if(!is.null(dnafile)){
                                               inputDNA(dnafile)
                                             }
                                           },
                                         
                                         # Method for inputting DNA sequences...
                                        inputDNA =
                                          function(input, format=NULL){
                                            DNA$InputDNA(input, format)
                                            userBlocks$initializePairsFromDNA(DNA)
                                            comparrisonSettings <<- ComparrisonSettings$new(DNA)
                                            if(triplets$tripletsGenerated()){
                                              warning("Loading a new sequence file into HybRIDS object. Deleting triplets and data from previous sequence file.")
                                              triplets <<- Triplets$new()
                                            }
                                            triplets$generateTriplets(DNA, comparrisonSettings, filesDirectory)
                                          },
                                         
                                         # Method for setting any parameter for any stage.
                                         setParameters =
                                           function(Step = NULL, ...){
                                             if(!any(Step == c("TripletGeneration", "SSAnalysis", "BlockDetection", "BlockDating", "Plotting"))){
                                               stop("You need to specify a valid analysis 'Step' to alter the paramerters of.\nThe steps are TripletGeneration, SSAnalysis, BlockDetection, BlockDating, and Plotting.")
                                             }
                                             Parameters <- list(...)
                                             if(Step == "TripletGeneration"){
                                               comparrisonSettings$setSettings(DNA, ...)
                                               triplets$generateTriplets(DNA, comparrisonSettings, filesDirectory)
                                             }
                                             if(Step == "SSAnalysis"){
                                               ssAnalysisSettings$setSettings(...) 
                                             }
                                             if(Step == "BlockDetection"){
                                               blockDetectionSettings$setSettings(...)
                                             }
                                             if(Step == "BlockDating"){
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
                                             if(Step == "Plotting"){
                                               plottingSettings$setSettings(...)
                                             }
                                           },
                                         
                                         # Method for analyzing the sequence similarity of triplets of sequences.
                                         analyzeSS = 
                                           function(tripletSelections = "NOT.SCANNED", replaceSettings = FALSE, ...){
                                             DNA$enforceDNA()
                                             if(length(list(...)) > 0){
                                               if(replaceSettings){
                                                 ssAnalysisSettings$setSettings(...)
                                                 settings <- ssAnalysisSettings
                                               } else {
                                                 settings <- ssAnalysisSettings$copy()
                                                 settings$setSettings(...)
                                               }
                                             } else {
                                               settings <- ssAnalysisSettings
                                             } 
                                             triplets$scanTriplets(tripletSelections, DNA, settings)
                                             message("Finished Sequence Similarity Analysis.")
                                           },
                                         
                                         # Method to execute the putative block finds.
                                         findBlocks =
                                           function(tripletSelections = "NOT.SEARCHED", replaceSettings = FALSE, ...){
                                             if(length(list(...)) > 0){
                                               if(replaceSettings){
                                                 blockDetectionSettings$setSettings(...)
                                                 settings <- blockDetectionSettings
                                               } else {
                                                 settings <- blockDetectionSettings$copy()
                                                 settings$setSettings(...)
                                               }
                                             } else {
                                               settings <- blockDetectionSettings
                                             }
                                             return(triplets$findBlocks(tripletSelections, settings))
                                             message("Finished potential blocks for all triplet selections.")
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
                                             settings <- plottingSettings
                                             if(length(list(...)) > 0){
                                               if(ReplaceParams){
                                                settings$setSettings(...) 
                                               } else {
                                                settings <- plottingSettings$copy()
                                                settings$setSettings(...)
                                               }
                                             }
                                             return(triplets$plotTriplets(Selections, plottingSettings))
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
                                         
                                         show = function() {
                                           cat("HybRIDS object:\n\n")
                                           DNA$show()
                                           cat("\n\n")
                                           comparrisonSettings$show()
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

