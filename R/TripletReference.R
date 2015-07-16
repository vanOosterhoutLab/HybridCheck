#' Reference class to store the results from sequence similarity analyses, and block detection runs for a given triplet.
#' @name SimilarityScan
#' @field TableFile A length 1 character vector storing the temporary filepath of the datatable of sequence similarity table.
#' @field Table Accessor function for the datatable of sequence similarity that is stored on temporary files.
#' @field WindowSizeUsed Single integer value, stores the size of the sliding window used for the scan.
#' @field StepSizeUsed Single integer value, stores the size of the sliding window used for the scan.
SimilarityScan <- setRefClass("SimilarityScan",
                              
                              fields = list(
                                TableFile = "character",
                                Table = function(value){
                                  if(missing(value)){
                                    read.table(TableFile)
                                    } else {
                                      write.table(value, file = TableFile)
                                    }
                                  },
                                WindowSizeUsed = "numeric",
                                StepSizeUsed = "numeric"),
                              
                              methods = list(
                                initialize =
                                  function(HCDir){
                                    "Method initializes the object, generates temporary filenames for the sequence similarity table."
                                    TableFile <<- tempfile(pattern = "SSTable", tmpdir = HCDir)
                                    blankTable()
                                  },
                                
                                blankTable =
                                  function(){
                                    "Method clears the SS analysis results table."
                                    Table <<- data.frame(WindowCenter = NA, WindowStart = NA, WindowEnd = NA,
                                                         ActualCenter = NA, ActualStart = NA, ActualEnd = NA,
                                                         AB = NA, AC = NA, BC = NA)
                                  },
                                
                                tableIsBlank =
                                  function(){
                                    "Returns TRUE, if the SS analysis table is blank and no results are contained in it."
                                    return(all(is.na(Table)))
                                  },
                                
                                finalize =
                                  function(){
                                    "Called when the object is destroyed, makes sure to delete the file saved in the system's temporary directory."
                                    unlink(TableFile)
                                  }
                                )
                              )

#' Reference class to store and manage triplet data.
#' @name Triplet
#' @field ContigNames character vector of length 3, stores the dna names of which the triplet was made.
#' @field ContigNumbers integer vector of length 3, stores the indecies of the sequences that made the triplet.
#'  
Triplet <- setRefClass("Triplet",
                       
                       fields = list(
                         ContigNames = "character",
                         ContigPairs = "list",
                         InformativeDNALength = "numeric",
                         FullDNALength = "numeric",
                         ScanData = "ANY",
                         Blocks = "list"
                         ),
                       
                       methods = list(
                         initialize = 
                           function(sequencenames, fullseqlength, HCDir){
                             "Initializer function, creates an instance of triplet."
                             ContigNames <<- sequencenames
                             FullDNALength <<- fullseqlength
                             ContigPairs <<- combn(ContigNames, 2, simplify=F)
                             ScanData <<- SimilarityScan$new(HCDir)
                           },
                         
                         readSettings =
                           function(dna, settings){
                             InformativeDNALength <<- unique(width(dna))
                             ScanData$StepSizeUsed <<- settings$StepSize
                             ScanData$WindowSizeUsed <<- settings$WindowSize
                           },
                         
                         noScanPerformed =
                           function(){
                             "Returns TRUE, if the SS analysis table is blank and the informative sites are not known. This is indicative that a scan of the file has not been done yet."
                             return(ScanData$tableIsBlank() && length(InformativeDNALength) == 0)
                           },
                         
                         blocksNotFound =
                           function(){
                             "Returns TRUE, if the tables for block are blank and no block findng has been done."
                             return(length(Blocks) == 0)
                           },
                         
                         blocksNotDated =
                           function(){
                             "Returns TRUE, if the blocks detected have not been tested for significance or had a divergence time estimated."
                             bools <- unlist(lapply(Blocks, function(y) all(unlist(lapply(y, function(x) all(is.na(x$SNPs)) && all(is.na(x$CorrectedSNPs)) && all(is.na(x$P_Value)) && all(is.na(x$P_Threshold)) && all(is.na(x$fiveAge)) && all(is.na(x$fiftyAge)) && all(is.na(x$ninetyFiveAge))   )))))
                             return(all(bools))
                           },
                         
                         # Method for putative block detection.
                         putativeBlockFind = 
                           function(parameters){
                             "DOCSTRING TO BE COMPLETE"
                             if(noScanPerformed()){stop("No sequence similarity scan data is available for this triplet - can't identify blocks.")}
                               if(parameters$AutoThresholds == TRUE) {
                                 message("Using the autodetect thresholds method...")
                                 message("Deciding on suitable thresholds...")
                                 thresholds <- autodetect.thresholds(ScanData, parameters)
                                 # Results in a list of thresholds for AB, AC and BC.
                               } else {
                                 thresholds <- list(parameters$ManualThresholds, parameters$ManualThresholds, parameters$ManualThresholds)
                               }
                               names(thresholds) <- unlist(lapply(combn(ContigNames, 2, simplify=F), function(x) paste(x, collapse=":")))
                               message("Now beginning Block Search...")
                             Blocks <<- lapply(1:3, function(i) block.find(ScanData$Table[,c( 1:6, 6+i )], thresholds[[i]]))
                             names(Blocks) <<- names(thresholds)
                           },
                         
                         blockDate =
                           function(dnaobj, parameters){
                             "Block Dating method, estimates the ages of blocks detected based on how many mutations are observed in a block and ."
                             message("Now dating blocks")
                             ab.blocks <- lapply(Blocks[[1]], function(x) date.blocks(x, dnaobj, parameters$MutationRate, ContigPairs[[1]], parameters$PValue, parameters$BonfCorrection, parameters$DateAnyway, parameters$MutationCorrection))
                             ac.blocks <- lapply(Blocks[[2]], function(x) date.blocks(x, dnaobj, parameters$MutationRate, ContigPairs[[2]], parameters$PValue, parameters$BonfCorrection, parameters$DateAnyway, parameters$MutationCorrection))
                             bc.blocks <- lapply(Blocks[[3]], function(x) date.blocks(x, dnaobj, parameters$MutationRate, ContigPairs[[3]], parameters$PValue, parameters$BonfCorrection, parameters$DateAnyway, parameters$MutationCorrection))
                             out.blocks <- list(ab.blocks, ac.blocks, bc.blocks)
                             names(out.blocks) <- names(Blocks)
                             Blocks <<- out.blocks
                           },
                         
                         tabulateBlocks = function(){
                           "Tabulates the blocks for a triplet."
                             message(paste("Tabulating blocks for the triplet", paste(ContigNames, collapse=":")))
                             # Check that the tables are present, if they aren't, turn them into blank data.frames.
                             temps <- lapply(1:3, function(i) do.call(rbind, Blocks[[i]]))
                             SS <- lapply(1:3, function(i) floor(as.numeric(rownames(temps[[i]]))))
                             pair <- lapply(1:3, function(i) rep(names(Blocks)[[i]], nrow(temps[[i]])))
                             temp2 <- do.call(rbind, temps)
                             if(blocksNotDated()){
                               temp2 <- cbind(temp2, data.frame(fiveAge = rep(NA, times=nrow(temp2)), fiftyAge = rep(NA, times=nrow(temp2)), ninetyfiveAge = rep(NA, times=nrow(temp2)), SNPnum = rep(NA, times=nrow(temp2)), PValue = rep(NA, times=nrow(temp2)), PThresh = rep(NA, times=nrow(temp2)), MeanAge = rep(NA, times=nrow(temp2)), CorrectedSNPs = rep(NA, times=nrow(temp2))))
                             }
                           temp2["SequencePair"] <- unlist(pair)
                           temp2["SequenceSimilarityThreshold"] <- unlist(SS)
                             if(nrow(temp2) > 0){
                               temp2["Triplet"] <- paste(ContigNames, collapse=":")
                             } else {
                               temp2["Triplet"] <- character()
                             }
                             return(temp2)
                         },
                         
                         plotTriplet = function(plottingSettings){
                           if("Lines" %in% plottingSettings$What && "Bars" %in% plottingSettings$What){
                             plottedbars <- plotBars(plottingSettings)
                             plottedlines <- plotLines(plottingSettings)
                             together <- arrangeGrob(plottedbars, plottedlines, ncol=1)
                             class(together) <- c("HC_LB_Plot", class(together))
                             return(together)
                           } else {
                             if("Lines" %in% plottingSettings$What){
                               plottedlines <- plotLines(plottingSettings)
                               class(plottedlines) <- c("HC_L_Plot", class(plottedlines))
                               return(plottedlines)
                             }
                             if("Bars" %in% plottingSettings$What){
                               plottedbars <- plotBars(plottingSettings)
                               class(plottedbars) <- c("HC_B_Plot", class(plottedbars))
                               return(plottedbars)
                             }
                           }
                         },
                         
                         plotLines =
                           function(plottingSettings){
                             "Method plots a lineplot using ggplot2 of the sequence similarity data from the scan."
                             if(noScanPerformed()){stop("No sequence similarity scan has been performed for this triplet.")}
                             combo <- unlist(lapply(combn(ContigNames, 2, simplify=FALSE), function(x) paste(x, collapse=":")))
                             data <- ScanData$Table
                             plotting.frame <- data.frame(basepos = rep(data$ActualCenter,3),
                                                           yvalues = c(data$AB, data$AC, data$BC),
                                                           factors = rep(1:3, each = nrow(data)))
                             plot <- ggplot(plotting.frame, aes(x=basepos, y=yvalues)) + geom_line(aes(colour=factor(factors)), show_guide=plottingSettings$Legends, size=0.8) +
                               ylim(0,100) + 
                               scale_colour_manual(name = "Pairwise Comparrisons", labels=c(combo[1], combo[2], combo[3]),values=c("yellow","purple","cyan")) +
                               xlab("Base Position") +
                               ylab("% Sequence Similarity")
                             plot <- applyPlottingParams(plot, plottingSettings, title = paste("Sequence Similarity Between Sequences for Triplet ", ContigNames[1], ":", ContigNames[2], ":", ContigNames[3], sep=""))
                             return(plot)
                           },
                         
                         plotBars =
                           function(plottingSettings){
                             "Method plots the heatmap based graphic of bars, from the sequence similarity scan data."
                             if(noScanPerformed()){stop("No sequence similarity scan has been performed for this triplet.")}
                             # Generate the reference colour palette.
                             colourPalette <- expand.grid(A = seq(0, 100, by = 1), B = seq(0, 100, by = 1))
                             colourPalette$RefA <- rgb(green = colourPalette$A, red = 100, blue = colourPalette$B, maxColorValue = 100)
                             colourPalette$RefB <- rgb(green = 100, red = colourPalette$A, blue = colourPalette$B, maxColorValue = 100)
                             colourPalette$RefC <- rgb(green = colourPalette$B, red = colourPalette$A, blue = 100, maxColorValue = 100)
                             # Now figure out the scale and data to go into each vertical bar: TODO - Put this in a function.
                             div <- FullDNALength / plottingSettings$MosaicScale
                             frame <- data.frame(bpstart = seq(from = 1, to = FullDNALength, by = div),
                                                 bpend = seq(from=div, to = FullDNALength, by = div)) 
                             frame$bpX <- round(frame$bpstart +  (div / 2))
                             scanTable <- ScanData$Table
                             AB <- round(apply(frame, 1, function(x) vertbar_create(scanTable, x, 7)))
                             AC <- round(apply(frame, 1, function(x) vertbar_create(scanTable, x, 8)))
                             BC <- round(apply(frame, 1, function(x) vertbar_create(scanTable, x, 9)))
                             rm(scanTable)
                             frame$AB <- AB
                             frame$AC <- AC
                             frame$BC <- BC
                             frame$X <- 1:nrow(frame)
                             rm(AB, AC, BC)
                             if(any(is.nan( as.numeric(frame$AB))) || any(is.nan( as.numeric(frame$AC))) || any(is.nan( as.numeric(frame$BC)))){
                               warning("\nNot a numbers (NaNs)! have been detected in the plotting frame.\n
The common cause of this is a small alignment or few informative sites in the data, 
with a too high MosaicScale parameter.\nThis usually happens at either end of the 
bars and the NaNs will be dealt with my filling them in black.\n\nTo get rid of them use a lower MosaicScale parameter.")
                             }
                             A_mix <- apply(frame, 1, function(x) col_deter(x[c(4,5)], colourPalette[,c(1,2,3)]))
                             B_mix <- apply(frame, 1, function(x) col_deter(x[c(4,6)], colourPalette[,c(1,2,4)]))
                             C_mix <- apply(frame, 1, function(x) col_deter(x[c(5,6)], colourPalette[,c(1,2,5)]))
                             frame$A_mix <- A_mix
                             frame$B_mix <- B_mix
                             frame$C_mix <- C_mix
                             rm(A_mix, B_mix, C_mix)
                             plottingFrame <- data.frame(X = frame$X, Y = rep(c(3, 2, 1), each = plottingSettings$MosaicScale), colour = c(frame$A_mix, frame$B_mix, frame$C_mix))
                             bars <- ggplot(plottingFrame, aes(x = X, y = as.factor(Y))) +
                               geom_raster(aes(fill = colour)) + scale_fill_identity() +
                               xlab("Approximate Base Position") +
                               ylab("Sequence Name") +
                               scale_x_continuous(breaks = c(seq(from = 1, to = plottingSettings$MosaicScale, by = plottingSettings$MosaicScale / 10), plottingSettings$MosaicScale), labels = c(frame$bpX[seq(from = 1, to = plottingSettings$MosaicScale, by = plottingSettings$MosaicScale / 10)], max(frame$bpX))) + 
                               scale_y_discrete(labels = c(ContigNames[3], ContigNames[2], ContigNames[1]))
                             bars <- applyPlottingParams(bars, plottingSettings, title = paste("Sequence Similarity Between Sequences for Triplet ", ContigNames[1], ":", ContigNames[2], ":", ContigNames[3], sep=""))
                             
                             if(plottingSettings$Legends == T){
                               legend <- readPNG(system.file("extdata/rgblegend.png", package="HybridCheck"), TRUE)
                               if (names(dev.cur()) == "windows"){
                                 # windows device doesn’t support semi-transparency so we’ll need
                                 # to flatten the image
                                 transparent <- legend[,,4] == 0
                                 legend <- as.raster(legend[,,1:3])
                                 legend[transparent] <- NA
                               }
                               legendgrob <- grid::rasterGrob(image=legend)
                               bars <- arrangeGrob(bars, legendgrob, widths = c(1, 0.13), ncol = 2)
                             }
                             return(bars)
                           }
                         )
                       )

vertbar_create <- function(sequenceSimilarityTable, plottingFrameRow, whichComparrison){
  bool1 <- sequenceSimilarityTable$ActualStart <= plottingFrameRow[2]
  bool2 <- plottingFrameRow[1] <= sequenceSimilarityTable$ActualEnd
  index <- which(bool1 == bool2)
  return(mean(sequenceSimilarityTable[index, whichComparrison]))
}

col_deter <- function(invalues, reference){
  if(any(is.nan(as.numeric(invalues)))){
    cols <- "#000000"
  } else {
    cols <- reference[(reference[, 1] == as.numeric(invalues[1])) & reference[,2] == as.numeric(invalues[2]), 3]
  }
  return(cols)
}




#' Reference class storing all triplets in a HC analysis.
#' @name Triplets
#' @description The Triplets reference class stores and manages operations over many Triplet objects.
Triplets <- setRefClass("Triplets",
                        
                        fields = list(
                          triplets = "list"
                          ),
                        
                        methods = list(
                          initialize =
                            function(dna){
                              "Initializes the object."
                              triplets <<- list()
                            },
                          
                          tripletsGenerated =
                            function(){
                              "Returns TRUE if the triplet generations have been prepared for."
                              return(length(triplets) > 0)
                            },
                          
                          matchNames =
                            function(selection){
                              if(!tripletsGenerated()){stop("No triplets have been prepared yet.")}
                              return(unlist(lapply(triplets, function(x) sum(selection %in% x$ContigNames))))
                            },
                          
                          deleteAllTriplets =
                            function(){
                              "Removes all current triplet data completely."
                              message("Deleting all triplets data.")
                              triplets <<- list()
                            },
                          
                          generateTriplets =
                            function(dna, csettings, basefile){
                              "Initializes all triplet objects based on the combination settings"
                              if(csettings$Modified){
                                if(tripletsGenerated()){
                                  deleteAllTriplets()
                                }
                                message("Initializing new triplets data.")
                                seqlength <- dna$getFullLength()
                                seqnames <- dna$getSequenceNames()
                                triplets <<- lapply(csettings$AcceptedCombinations, function(x) Triplet$new(sequencenames = c(x[1], x[2], x[3]), fullseqlength = seqlength, basefile))
                                csettings$Modified <- FALSE
                              }
                            },
                          
#                           updateTriplets =
#                             function(dna, settings, basefile){
#                               got <- getAllNames()
#                               to.remove <- which(!unlist(lapply(got, function(x) x %in% settings$AcceptedCombinations)))
#                               to.add <- which(!unlist(lapply(settings$AcceptedCombinations, function(x) x %in% got)))
#                               triplets <<- triplets[-to.remove]
#                               triplets <<- append(triplets, makeTriplets(settings$AcceptedCombinations[to.add], dna, basefile))
#                             },
#                           
#                           makeTriplets =
#                             function(selection, dna, basefile){
#                               seqlength <- dna$getFullLength()
#                               return(lapply(selection, function(x) Triplet$new(sequencenumbers = x, sequences = c(x[1], x[2], x[3]), fullseqlength = seqlength, basefile)))
#                             },
                          
                          scanTriplets =
                            function(tripletSelections, dna, scansettings){
                              if(!tripletsGenerated()){stop("No triplets have been prepared yet.")}
                              tripletsToScan <- getTriplets(tripletSelections)
                              for(tripletToScan in tripletsToScan){
                                scan.similarity(dna, tripletToScan, scansettings)
                              }
                            },

                          findBlocks =
                            function(tripletSelections, findSettings){
                              if(!tripletsGenerated()){stop("No triplets have been prepared yet.")}
                              tripletsToFindIn <- getTriplets(tripletSelections)
                              for(triplet in tripletsToFindIn){
                                triplet$putativeBlockFind(findSettings)
                              }
                            },

                          dateBlocks =
                            function(tripletSelections, dateSettings, dna){
                              if(!tripletsGenerated()){stop("No triplets have been prepared yet.")}
                              tripletsToDate <- getTriplets(tripletSelections)
                              for(triplet in tripletsToDate){
                                triplet$blockDate(dna, dateSettings)
                              }
                            },
                          
                          plotTriplets = function(tripletSelections, plotSettings){
                            tripletsToPlot <- getTriplets(tripletSelections)
                            return(lapply(tripletsToPlot, function(x) x$plotTriplet(plotSettings)))
                          },

                          tabulateBlocks = function(tripletSelections, neat){
                            if(!tripletsGenerated()){stop("No triplets have been prepared yet.")}
                            tripletsToTabulate <- getTriplets(tripletSelections)
                            listedTabulates <- lapply(tripletsToTabulate, function(x) x$tabulateBlocks())
                            collected <- do.call(rbind, listedTabulates)
                            output <- data.frame(collected$Triplet, collected$SequencePair, collected$SequenceSimilarityThreshold, collected$Length,
                                                 collected$First, collected$Last, collected$FirstBP, collected$LastBP, collected$ApproxBpLength, collected$SNPs, collected$CorrectedSNPs, collected$fiveAge, collected$fiftyAge,
                                                 collected$ninetyFiveAge, collected$P_Value, collected$P_Threshold)
                            if(neat){
                              output <- output[,-c(4,5,6)]
                              names(output) <- c("Triplet", "Sequence_Pair","Sequence_Similarity_Threshold","First_BP_Position", "Last_BP_Position", "Approximate_Length_BP", "Number_of_SNPs", "Corrected_Number_of_SNPs", "p=0.05_Age", "p=0.5_Age","p=0.95_Age","P_Value", "P_Thresh")
                            }
                            return(output)
                          },
                          
                          show =
                            function(){
                              cat("A total of ")
                            },

                          tripletCombinations =
                            function(){
                              return(unlist(lapply(triplets, function(x){paste0(x$ContigNames, collapse=", ")})))
                            },
                          
                          htmlSummary =
                            function(){
                              "Prints a HTML summary of all the different triplet combinations."
                              tripnames <- paste0(lapply(triplets, function(x){paste0(paste0(x$ContigNames, collapse=", "), " <br> ", collapse="")}), collapse="")
                              output <- paste0("<h2>TripletCombinations</h2>")
                              return(tripnames)
                            },
                          
                          getTriplets =
                            function(selections){
                              "Returns a list of references to triplets according to user selection."
                              if(!is.list(selections)){
                                selections <- list(selections)
                              }
                              selections <- unique(selections)
                              if(!is.null(selections) && length(selections) > 0){
                                ind <- numeric()
                                if("ALL" %in% selections){
                                  ind <- 1:length(triplets)
                                } else {
                                  if("NOT.SCANNED" %in% selections){
                                    ind <- c(ind, which(unlist(lapply(triplets, function(x) x$noScanPerformed()))))
                                    selections <- selections[which(selections != "NOT.SCANNED")]
                                  }
                                  if("NOT.SEARCHED" %in% selections){
                                    ind <- c(ind, which(unlist(lapply(triplets, function(x) x$blocksNotFound()))))
                                    selections <- selections[which(selections != "NOT.SEARCHED")]
                                  }
                                  if("NOT.DATED" %in% selections){
                                    ind <- c(ind, which(unlist(lapply(triplets, function(x) x$blocksNotDated()))))
                                    selections <- selections[which(selections != "NOT.DATED")]
                                  }
                                  if(any(unlist(lapply(selections, length)) != 3)){stop("Selections must provide a vector of 3 sequence names.")}
                                  if(any(unlist(lapply(selections, function(x) !is.character(x))))){stop("Selections must be of class character.")}
                                  allNames <- do.call(rbind, getAllNames())
                                  ind <- c(ind, unlist(lapply(selections, function(x) which(allNames[,1] %in% x & allNames[,2] %in% x & allNames[,3] %in% x))))
                                }
                                ind <- unique(ind)
                                return(triplets[ind])
                              } else {
                                stop("No selection of triplet was provided.")
                              }
                            },
                          
                          getAllNames = 
                            function(){
                              "Returns the names of the sequences in each triplet as a list."
                              return(lapply(triplets, function(x) x$ContigNames))
                            },
                          
                          getAllIndexes =
                            function(){
                              "Returns the indexes of the sequences (according to their rows in the HC sequence object) in each triplet as a list."
                              return(lapply(triplets, function(x) x$ContigIndexes))
                            }
                          )
                        )
