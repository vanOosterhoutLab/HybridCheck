#' @docType package
#' ...
#' @useDynLib HybRIDS
#' @import ape ggplot2 grid gridExtra png
NULL

#' A Reference Class for managing a HybRIDS analysis.
#' @name HybRIDS
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
                          blockDatingSettings = "ANY",
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
                                             #BlockDatingParams <<- list(MutationRate = 10e-08, PValue = 0.005, BonfCorrection = TRUE, DateAnyway = FALSE)
                                             blockDatingSettings <<- BlockDatingSettings$new()
                                             
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
                                         
                                        showParameters =
                                          function(Step = NULL){
                                            for(i in Step){
                                              if(i == "TripletGeneration"){
                                                comparrisonSettings$show()
                                              }
                                              if(i == "SSAnalysis"){
                                                ssAnalysisSettings$show()
                                              }
                                              if(i == "BlockDetection"){
                                                blockDetectionSettings$show()
                                              }
                                              if(i == "BlockDetection"){
                                                blockDatingSettings$show()
                                              }
                                              if(i == "Plotting"){
                                                plottingSettings$show()
                                              }
                                              cat('\n\n')
                                            }
                                          },
                                        
                                         # Method for setting any parameter for any stage.
                                         setParameters =
                                           function(Step = NULL, ...){
                                             if(!any(Step == c("TripletGeneration", "SSAnalysis", "BlockDetection", "BlockDating", "Plotting")) || length(Step) != 1){
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
                                               blockDatingSettings$setSettings(...)
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
                                             triplets$findBlocks(tripletSelections, settings)
                                             message("Finished finding potential blocks for all triplet selections.")
                                           },
                                         
                                         # Method to Date the blocks found.
                                         dateBlocks =
                                           function(tripletSelections = "NOT.DATED", replaceSettings = FALSE, ...){
                                             if(length(list(...)) > 0){
                                               if(replaceSettings){
                                                 blockDatingSettings$setSettings(...)
                                                 settings <- blockDatingSettings
                                               } else {
                                                 settings <- blockDatingSettings$copy()
                                                 settings$setSettings(...)
                                               }
                                             } else {
                                               settings <- blockDatingSettings
                                             }
                                             triplets$dateBlocks(tripletSelections, settings, DNA)
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
                                           function(Selection = "ALL", Neat = TRUE) {
                                             output <- triplets$tabulateBlocks(Selection, Neat)
                                             class(output) <- c(class(output), "HybRIDStable")   
                                             return(output)
                                           },
                                         
                                         show = function() {
                                           cat("HybRIDS object:\n\n")
                                           DNA$show()
                                           cat("\n\n\n")
                                           comparrisonSettings$show()
                                           cat("\n\n\n")
                                           ssAnalysisSettings$show()
                                           cat("\n\n\n")
                                           blockDetectionSettings$show()
                                           cat("\n\n\n")
                                           blockDatingSettings$show()
                                           cat("\n\nPlotting Settings are not shown because of the number of settings.")
                                           cat("\nuse showParameters('Plotting') to view them.")
                                         },
                                         
                                         addUserBlock = function(selection, firstbp, lastbp){
                                           userBlocks$addBlock(firstbp, lastbp, selection)
                                         },
                                        
                                         clearUserBlocks = function(selection){
                                           userBlocks$blankBlocks(selection)
                                         },
                                        
                                         dateUserBlocks = function(){
                                           userBlocks$dateBlocks(DNA, blockDatingSettings)
                                         },
                                        
                                        tabulateUserBlocks = function(){
                                          return(userBlocks$tabulateBlocks())
                                        }
                          )
                        )

