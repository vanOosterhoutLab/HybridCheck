#' @docType package
#' ...
#' @useDynLib HybRIDS
#' @import methods ape ggplot2 gridExtra Biostrings
#' @importFrom png readPNG
NULL

#' A Reference Class for managing a HybRIDS analysis.
#' @name HybRIDS
#' @export HybRIDS
#' @exportClass HybRIDS
#' @field DNA A HybRIDSdna reference object.
#' @field FTTmodule A FTTester object, controls and stores results of the Four Taxon Tests.
#' @field comparrisonSettings A ComparrisonSettings object.
#' @field ssAnalysisSettings A SSAnalysisSettings object.
#' @field blockDetectionSettings A BlockDetectionSettings object.
#' @field blockDatingSettings A BlockDatingSettings object.
#' @field plottingSettings A PlottingSettings object.
#' @field triplets A Triplets object.
#' @field userBlocks A UserBlocks object.
#' @field filesDirectory Character - the root directory where all temporary files 
#' used by this object are found.
HybRIDS <- setRefClass("HybRIDS",
                        
                        fields = list( 
                          DNA = "ANY",
                          FTTmodule = "ANY",
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
                                             "Create HybRIDS object with default values for fields.
                                             The path to the FASTA file can be provided."
                                             filesDirectory <<- tempdir()
                                             DNA <<- HybRIDSseq$new()
                                             FTTmodule <<- FTTester$new()
                                             ssAnalysisSettings <<- SSAnalysisSettings$new()
                                             blockDetectionSettings <<- BlockDetectionSettings$new()
                                             blockDatingSettings <<- BlockDatingSettings$new()
                                             plottingSettings <<- PlottingSettings$new()
                                             userBlocks <<- UserBlocks$new()
                                             triplets <<- Triplets$new()
                                             if(!is.null(dnafile)){
                                               inputDNA(dnafile)
                                             }
                                           },
                                         
                                         # Method for inputting DNA sequences...
                                        inputDNA =
                                          function(input){
                                            "Loads FASTA file or Biostrings object into the session."
                                            DNA$InputDNA(input)
                                            userBlocks$initializePairsFromDNA(DNA)
                                            comparrisonSettings <<- ComparrisonSettings$new(DNA, FTTmodule)
                                            if(triplets$tripletsGenerated()){
                                              warning("Loading a new sequence file into HybRIDS object. Deleting triplets and data from previous sequence file.")
                                              triplets <<- Triplets$new()
                                            }
                                            setPopulations()
                                            triplets$generateTriplets(DNA, comparrisonSettings, filesDirectory)
                                          },
                                        
                                        sequencesLoaded =
                                          function(){
                                            "Returns TRUE if dna sequences have been loaded."
                                            return(DNA$hasDNA())
                                          },
                                        
                                        sequenceNames =
                                          function(){
                                            "Returns a character vector of all the sequence names loaded."
                                            DNA$getSequenceNames()
                                          },
                                        
                                        setPopulations =
                                          function(pops = NULL){
                                            "Sets up population definitions. If pops is NULL,
                                            each sequence loaded is assumed to be from one population."
                                            if(is.null(pops)){
                                              DNA$oneSeqOnePop()
                                            } else {
                                              DNA$setPopulations(pops)
                                            }
                                          },
                                        
                                        prepareFourTaxonTests =
                                          function(taxa = NULL){
                                            "Prepares space for four taxon tests."
                                            if(is.null(taxa)){
                                              FTTmodule$autoTaxaCombos(DNA)
                                            } else {
                                              FTTmodule$manualTaxaCombos(taxa, DNA)
                                            }
                                            FTTmodule$generateFTTs(filesDirectory)
                                          },
                                        
                                        runFourTaxonTests =
                                          function(selections = "NOT.TESTED", numberOfBlocks = NULL, blockLength = NULL){
                                            "Runs the four taxon tests that have been prepared.
                                            You must provide either a number of blocks to use for the jackknife, 
                                            or the size of the blocks to use for the jackknife, 
                                            the exact size and number of blocks to use is then automatically worked out.
                                            selections can be 'NOT.TESTED', 'ALL', 'TESTED', or a list of character vectors of length 4, each 
                                            denoting a four taxon test by the names of the populations involved."
                                            FTTmodule$runFTTests(selections, DNA, numberOfBlocks, blockLength)
                                            if(1L %in% comparrisonSettings$Method){
                                              comparrisonSettings$decideAcceptedTriplets(DNA, FTTmodule)
                                              triplets$generateTriplets(DNA, comparrisonSettings, filesDirectory)
                                            }
                                          },
                                        
                                        tabulateFourTaxonTests =
                                          function(selections = "ALL", neat = TRUE, global = TRUE){
                                            "Fetches the results of the four taxon tests, returns a dataframe.
                                            selections can be 'ALL', 'TESTED', or a list of character vectors of length 4, each 
                                            denoting a four taxon test by the names of the populations involved."
                                            FTTmodule$getResults(selections, neat = neat)
                                          },
                                         
                                        showParameters =
                                          function(Step = NULL){
                                            "Displays to the R console, the settings of each HybRIDS analysis stage."
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
                                        
                                         setParameters =
                                           function(Step = NULL, ...){
                                             "Set parameters for a single stage of the analysis.
                                             Specify a single stage for Step.
                                             The steps are: 'FTT', 'TripletGeneration',
                                             'SSAnalysis', 'BlockDetection', 'BlockDating',
                                             'Plotting'."
                                             if(!any(Step == c("TripletGeneration", "SSAnalysis", "BlockDetection", "BlockDating", "Plotting")) || length(Step) != 1){
                                               stop("You need to specify a valid analysis 'Step' to alter the paramerters of.\nThe steps are TripletGeneration, SSAnalysis, BlockDetection, BlockDating, and Plotting.")
                                             }
                                             Parameters <- list(...)
                                             if(Step == "TripletGeneration"){
                                               comparrisonSettings$setSettings(DNA, FTTmodule, ...)
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
                                             if(Step == "FTT"){
                                               FTTmodule$setSettings(...)
                                               FTTmodule$generateFTTs(filesDirectory)
                                             }
                                           },
                                         
                                         analyzeSS = 
                                           function(tripletSelections = "NOT.SCANNED", replaceSettings = FALSE, ...){
                                             "Analyzing the sequence similarity of triplets of sequences.
                                             tripletSelections can be 'NOT.SCANNED', 'ALL', 'SCANNED',
                                             or a list of character vectors of length 3, each specifying a sequence triplet
                                             by the names of the sequences involved. If replaceSettings is TRUE, then
                                             optional settings passed for this analysis will apply to subsequent calls
                                             of this function."
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
                                           function(Selections = "ALL", ReplaceParams = TRUE, ...){
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

