PlottingSettings <- setRefClass("PlottingSettings",
                                
                                fields = list(
                                  What = "character",
                                  PlotTitle = "logical",
                                  CombinedTitle = "logical",
                                  TitleSize = "numeric",
                                  TitleFace = "character",
                                  TitleColour = "character",
                                  XLabels = "logical",
                                  YLabels = "logical",
                                  XTitle = "logical",
                                  XTitleFontSize = "numeric",
                                  XTitleColour = "character",
                                  XLabelSize = "numeric",
                                  XLabelColour = "character",
                                  YTitle = "logical",
                                  YTitleFontSize = "numeric",
                                  YTitleColour = "character",
                                  YLabelSize = "numeric",
                                  YLabelColour = "character",
                                  Legends = "logical",
                                  LegendFontSize = "numeric",
                                  MosaicScale = "numeric"),
                                
                                methods = list(
                                  initialize = function(){
                                    What <<- c("Bars", "Lines")
                                    PlotTitle <<- TRUE
                                    CombinedTitle <<- FALSE 
                                    TitleSize <<- 14
                                    TitleFace <<- "bold"
                                    TitleColour <<- "black"
                                    XLabels <<- TRUE
                                    YLabels <<- TRUE
                                    XTitle <<- TRUE
                                    XTitleFontSize <<- 12
                                    XTitleColour <<- "black"
                                    XLabelSize <<- 10
                                    XLabelColour <<- "black"
                                    YTitle <<- TRUE
                                    YTitleFontSize <<- 12
                                    YTitleColour <<- "black"
                                    YLabelSize <<- 10
                                    YLabelColour <<- "black"
                                    Legends <<- TRUE
                                    LegendFontSize <<- 12
                                    MosaicScale <<- 500
                                  },
                                  
                                  setWhat = function(value){
                                    if(!all(value == "Lines" | value == "Bars") || length(value) > 2){stop("Input value must be 'Bars' and/or 'Lines'")}
                                    What <<- value
                                  },
                                  
                                  textSummary =
                                    function(){
                                      return(paste('Settings for plotting recombination signal in triplets:\n',
                                                   '-------------------------------------------------------\n\n',
                                                   'What graphs will be plotted (What): ', paste(What, collapse=", "),
                                                   '\n\nWhether to include a title with the plot. (PlotTitle): ',
                                                   PlotTitle,
                                                   '\n\nOne overall title for a Lines and Bars plot on the same canvas (CombinedTitle): ',
                                                   CombinedTitle,
                                                   '\n\nFont size of titles (TitleSize): ', TitleSize,
                                                   '\n\nFace of the title (TitleFace): ', TitleFace,
                                                   '\n\nColour of the plot titles (TitleColour): ', TitleColour,
                                                   '\n\nInclude labels of the x axes in plots (XLabels): ', XLabels,
                                                   '\n\nSize of the x axes labels (XLabelSize): ', XLabelSize,
                                                   '\n\nColiur of the x axes labels (XLabelColour): ', XLabelColour,
                                                   '\n\nInclude the title of the x axis in plots (XTitle): ', XTitle,
                                                   '\n\nFont size of the x axis title (XTitleFontSize): ', XTitleFontSize,
                                                   '\n\nColour of the x axis title (XTitleColour): ', XTitleColour,
                                                   '\n\nInclude labels of the y axes in plots (YLabels): ', YLabels,
                                                   '\n\nSize of the y axes labels (YLabelSize): ', YLabelSize,
                                                   '\n\nColiur of the y axes labels (YLabelColour): ', YLabelColour,
                                                   '\n\nInclude the title of the y axis in plots (YTitle): ', YTitle,
                                                   '\n\nFont size of the y axis title (YTitleFontSize): ', YTitleFontSize,
                                                   '\n\nColour of the y axis title (YTitleColour): ', YTitleColour,
                                                   '\n\nInclude legends for plots (Legends): ', Legends,
                                                   '\n\nFont size for legends (LegendFontSize): ', LegendFontSize,
                                                   '\n\nMosaicScale: ', MosaicScale      
                                      ))
                                    },
                                  
                                  show = function(){
                                    cat(textSummary())
                                  },
                                  
                                  setSettings = function(...){
                                    settings <- list(...)
                                    parameters <- names(settings)
                                    for(i in 1:length(settings)){
                                      switch(
                                        parameters[i],
                                        What = setWhat(settings[[i]]),
                                        PlotTitle = PlotTitle <<- settings[[i]],
                                        CombinedTitle = CombinedTitle <<- settings[[i]],
                                        TitleSize = TitleSize <<- settings[[i]],
                                        TitleFace = TitleFace <<- settings[[i]],
                                        TitleColour = TitleColour <<- settings[[i]],
                                        XLabels = XLabels <<- settings[[i]],
                                        YLabels = YLabels <<- settings[[i]],
                                        XTitle = XTitle <<- settings[[i]],
                                        XTitleFontSize = XTitleFontSize <<- settings[[i]],
                                        XTitleColour = XTitleColour <<- settings[[i]],
                                        XLabelSize = XLabelSize <<- settings[[i]],
                                        XLabelColour = XLabelColour <<- settings[[i]],
                                        YTitle = YTitle <<- settings[[i]],
                                        YTitleFontSize = YTitleFontSize <<- settings[[i]],
                                        YTitleColour = YTitleColour <<- settings[[i]],
                                        YLabelSize = YLabelSize <<- settings[[i]],
                                        YLabelColour = YLabelColour <<- settings[[i]],
                                        Legends = Legends <<- settings[[i]],
                                        LegendFontSize = LegendFontSize <<- settings[[i]],
                                        MosaicScale = MosaicScale <<- settings[[i]]
                                      )
                                    }
                                  }
                                )
)

#' A reference class to represent settings for triplet generation.
#' @name ComparrisonSettings
#' @field Method An integer vector of length 1.
#' @field DistanceThreshold A numeric vector of length 1.
#' @field PartitionStrictness An integer vector of length 1.
#' @field Refine A logical vector of length 1.
ComparrisonSettings <- setRefClass("ComparrisonSettings",
                                   
                                   fields = list(
                                     SeqNames = "character",
                                     Method = "integer",
                                     DistanceThreshold = "numeric",
                                     PartitionStrictness = "integer",
                                     partiallySignificant = "logical",
                                     TripletCombinations = "list",
                                     AcceptedCombinations = "list",
                                     Modified = "logical"),
                                   
                                   methods = list(
                                     initialize = 
                                       function(dna, ftt){
                                         "Creates the object and sets all parameters to their default."
                                         Method <<- 1L
                                         DistanceThreshold <<- 0.01
                                         PartitionStrictness <<- 2L
                                         partiallySignificant <<- FALSE
                                         SeqNames <<- dna$getSequenceNames()
                                         TripletCombinations <<- combn(dna$getSequenceNames(), 3, simplify=FALSE)
                                         AcceptedCombinations <<- list()
                                         decideAcceptedTriplets(dna, ftt)
                                       },
                                     
                                     changeMethod =
                                       function(value){
                                         value <- unique(value)
                                         if(!is.integer(value)){stop("You must enter integer values between 1 and 4.")}
                                         if(any(value > 4L)){stop("4L is the maximum value allowed.")}
                                         if(any(value < 1L)){stop("1L is the minimum value allowed.")}
                                         if(3L %in% value && 4L %in% value){stop("Choose either method 3L or 4L, not both.")}
                                         Method <<- value
                                       },
                                     
                                     setDistanceThreshold =
                                       function(value){
                                         if(length(value) != 1 || !(value > 0 && value < 1)){stop("You must give a single double value, between 0 and 1, as a distance threshold.")}
                                         DistanceThreshold <<- value
                                       },
                                     
                                     setPartitionStrictness =
                                       function(value){
                                         if(length(value) != 1 || !any(value == c(1L, 2L))){stop("You must enter a single integer value of 1 or 2 as a PartitionStrictness.")}
                                         PartitionStrictness <<- value
                                       },
                                     
                                     setPartiallySignificant =
                                       function(value){
                                         if(length(value) > 1){stop("Only provide one logical value.")}
                                         partiallySignificant <<- value
                                       },
                                     
                                     decideAcceptedTriplets =
                                       function(dna, ftt){
                                         AcceptedCombinations <<- TripletCombinations
                                         rejects <- c()
                                         if(hasMultipleCombinations()){
                                           if(1L %in% Method){
                                             message(" - Deciding triplets based on results of Four Taxon Tests.")
                                             if(partiallySignificant){
                                               rejects <- c(rejects, fttDescision(ftt, "PART.SIGNIFICANT", TripletCombinations, dna))
                                             } else {
                                               message("\t- Using tests that are globally significant.")
                                               rejects <- c(rejects, fttDescision(ftt, "SIGNIFICANT", TripletCombinations, dna))
                                             }
                                           }
                                           if(2L %in% Method && dna$hasPopulations()){
                                             rejects <- c(rejects, groupDescision(dna$Populations, TripletCombinations, PartitionStrictness)) 
                                           }
                                           if(3L %in% Method || 4L %in% Method){                                                    
                                             rejects <- c(rejects, distanceDescision(dna, Method, DistanceThreshold, TripletCombinations))
                                           }
                                         } else {
                                           warning("There is only one comparrison possible - presumably only 3 sequences are present.")
                                         }
                                         if(!is.null(rejects) && length(rejects) > 0){
                                           AcceptedCombinations <<- AcceptedCombinations[-rejects]
                                         }
                                         Modified <<- TRUE
                                       },
                                     
                                     hasTripletCombinations =
                                       function(){
                                         return(length(TripletCombinations) > 0)
                                       },
                                     
                                     hasMultipleCombinations =
                                       function(){
                                         return(length(TripletCombinations) > 1)
                                       },
                                     
                                     numberOfTripletCombinations =
                                       function(){
                                         return(length(TripletCombinations))
                                       },
                                     
                                     numberOfAcceptedCombinations =
                                       function(){
                                         return(length(AcceptedCombinations))
                                       },
                                     
                                     textSummary =
                                       function(){
                                         "Creates a character vector of a summary of the comparrison settings."
                                         return(paste('Settings for Sequence Scan Combinations:\n',
                                                      '----------------------------------------\n',
                                                      'Triplet Generation Method (Method): ', paste0(Method, collapse = ", "),
                                                      '\n\nDistance Threshold for excluding triplets with\n\ttoo similar sequences (DistanceThreshold): ', DistanceThreshold,
                                                      '\n\nHow many sequences from the same partition are\n\tallowed in a triplet (PartitionStrictness): ', PartitionStrictness, 
                                                      '\n\nA total of ', numberOfAcceptedCombinations(), ' triplets will be compared.', 
                                                      sep=""))
                                       },
                                     
                                     printAcceptedCombinations =
                                       function(){
                                         return(lapply(AcceptedCombinations, function(x) paste0(x, collapse=", ")))
                                       },
                                     
                                     htmlCombinations =
                                       function(){
                                         return(paste0(lapply(AcceptedCombinations, function(x) paste0(x, collapse=", ")), collapse="<br>"))
                                       },
                                     
                                     show = function(){
                                       "Prints a summary of the object to console."
                                       cat(textSummary())
                                     },
                                     
                                     setSettings =
                                       function(dna, ftt, ...){
                                         settings <- list(...)
                                         parameters <- names(settings)
                                         for(i in 1:length(settings)){
                                           if(parameters[i] == "Method"){
                                             changeMethod(settings[[i]])
                                           }
                                           if(parameters[i] == "DistanceThreshold"){
                                             setDistanceThreshold(settings[[i]])
                                           }
                                           if(parameters[i] == "PartitionStrictness"){
                                             setPartitionStrictness(settings[[i]])
                                           }
                                           if(parameters[i] == "partiallySignificant"){
                                             setPartiallySignificant(settings[[i]])
                                           }
                                         }
                                         decideAcceptedTriplets(dna, ftt)
                                       }
                                   )
)

fttDescision <- function(ftt, significanceStatement, combos, dna){
  significantCombos <- ftt$getFTTs(significanceStatement)
  if(length(significantCombos) > 0){
    popnames <- lapply(significantCombos, function(x) x$getPops())
    matches <- sort(unique(unlist(lapply(popnames, function(x){
      seqnames <- dna$Populations[x]
      which(unlist(lapply(combos, function(y){
        return(sum(any(y %in% seqnames[[1]]), any(y %in% seqnames[[2]]),
                   any(y %in% seqnames[[3]]), any(y %in% seqnames[[4]])) == 3)
      })))
    }))))
    rejects <- which(!1:length(combos) %in% matches)
  } else {
    warning("No four taxon tests were found to be significant, either none were significant or no test has been performed.")
    rejects <- numeric()
  }
  return(rejects)
}

groupDescision <- function(populations, combos, thresh){
  message(" - Generating triplets to find recombination in sequences, between partitions.")
  matches <- unlist(lapply(populations, function(x){
    which(unlist(lapply(combos, function(y){
      length(which(x %in% y)) > thresh
    })))
  }))
  return(matches)
}

distanceDescision <- function(dna, method, thresh, combos){
  rejects <- c()
  distances <- dist.dna(as.DNAbin(dna$FullSequence), model = "raw")
  seqpairs <- combn(dna$getSequenceNames(), 2, simplify=FALSE)
  if(3L %in% method){
    rejectiondistances <- seqpairs[which(distances < thresh)]
  }
  if(4L %in% method){
    distances_density <- density(distances)
    Lows <- cbind(distances_density$x[which(diff(sign(diff(distances_density$y))) == 2)], distances_density$y[which(diff(sign(diff(distances_density$y))) == 2)])
    Lowest <- Lows[which(Lows[,1] == min(Lows[,1])),]
    rejectiondistances <- seqpairs[which(distances < Lowest[1])]
  }
  for(i in 1:length(combos)){
    for(n in 1:length(rejectiondistances)){
      if(all(rejectiondistances[[n]] %in% combos[[i]])){
        rejects <- c(rejects, i)
        break
      }
    }
  }
  return(rejects)
}

applyPlottingParams <- function(plot, parameters, title = ""){
  if(parameters$PlotTitle){
    plot <- plot + 
      ggtitle(title) +
      theme(title = element_text(size = parameters$TitleSize, colour = parameters$TitleColour, face = parameters$TitleFace),
            legend.text = element_text(size = parameters$LegendFontSize))
  } else {
    plot <- plot + theme(title = element_blank(),
                         legend.text = element_text(size = parameters$LegendFontSize))
  }
  if(parameters$XTitle){
    plot <- plot + theme(axis.title.x = element_text(size = parameters$XTitleFontSize, colour = parameters$XTitleColour))
  } else {
    plot <- plot + theme(axis.title.x = element_blank())
  }
  if(parameters$YTitle){
    plot <- plot + theme(axis.title.y = element_text(size = parameters$YTitleFontSize, colour = parameters$YTitleColour))
  } else {
    plot <- plot + theme(axis.title.y = element_blank())
  }
  if(parameters$XLabels){
    plot <- plot + theme(axis.text.x = element_text(size = parameters$XLabelSize, colour = parameters$XLabelColour))
  } else {
    plot <- plot + theme(axis.text.x = element_blank())
  }
  if(parameters$YLabels){
    plot <- plot + theme(axis.text.y = element_text(size = parameters$YLabelSize, colour = parameters$YLabelColour))
  } else {
    plot <- plot + theme(axis.text.y = element_blank())
  }
  return(plot)
}

SSAnalysisSettings <- setRefClass("SSAnalysisSettings",
                                  
                                  fields = list(
                                    WindowSize = "integer",
                                    StepSize = "integer"
                                  ),
                                  
                                  methods = list(
                                    initialize =
                                      function(winSize = NULL, stepSize = NULL){
                                        if(is.null(winSize)){
                                          WindowSize <<- 100L
                                        } else {
                                          WindowSize <<- winSize
                                        }
                                        if(is.null(stepSize)){
                                          StepSize <<- 1L
                                        } else {
                                          StepSize <<- stepSize
                                        }
                                      },
                                    
                                    getWindowSize =
                                      function(){
                                        return(WindowSize)
                                      },
                                    
                                    getStepSize =
                                      function(){
                                        return(StepSize)
                                      },
                                    
                                    setWindowSize =
                                      function(value){
                                        if(length(value) != 1 || !is.integer(value)){stop("Error: Provide one integer value as a Window Size.")}
                                        WindowSize <<- value
                                      },
                                    
                                    setStepSize =
                                      function(value){
                                        if(length(value) != 1 || !is.integer(value)){stop("Error: Provide one integer value as a Step Size.")}
                                        StepSize <<- value
                                      },
                                    
                                    setSettings =
                                      function(...){
                                        settings <- list(...)
                                        parameters <- names(settings)
                                        for(i in 1:length(settings)){
                                          if(parameters[i] == "WindowSize"){
                                            setWindowSize(settings[[i]])
                                          }
                                          if(parameters[i] == "StepSize"){
                                            setStepSize(settings[[i]])
                                          }
                                        }
                                      },
                                    
                                    textSummary =
                                      function(){
                                        return(paste('Settings for sliding window scan of recombination signal:\n',
                                                     '---------------------------------------------------------\n',
                                                     'Size of the sliding window in base pairs (WindowSize): ',
                                                     getWindowSize(),
                                                     '\n\nNumber of base pairs to move the sliding window on\n\teach iteration of the scan (StepSize): ', 
                                                     getStepSize(), sep=""))
                                      },
                                    
                                    show = 
                                      function(){
                                        cat(textSummary())
                                      }
                                  )
)


# Functions and classes for detecting blocks.
#' Reference type for storing the settings for block detection.
#' @name BlockDetectionSettings
#' @field ManualThresholds A numeric value between 1 and 100, the percentage raw sequence similarity a region has to reach, before it is identified as a block.
BlockDetectionSettings <- setRefClass("BlockDetectionSettings",
                                      
                                      fields = list(
                                        ManualThresholds = "numeric",
                                        AutoThresholds = "logical",
                                        ManualFallback = "logical",
                                        SDstringency = "numeric"
                                      ),
                                      
                                      methods = list(
                                        initialize = function(){
                                          "Initializes the settingsobject."
                                          ManualThresholds <<- 90
                                          AutoThresholds <<- TRUE
                                          ManualFallback <<- TRUE
                                          SDstringency <<- 2
                                        },
                                        
                                        setManualThresholds =
                                          function(values){
                                            "Checks input values for changing the manual sequence similarity thresholds for block detection and sets the parmeter."
                                            if(any(values > 100) || any(values < 0)){stop("Enter a numeric value between 1 and 100.")}
                                            ManualThresholds <<- values
                                          },
                                        
                                        setSDstringency =
                                          function(value){
                                            "Checks input values for changing the SD stringency parameter for block detection and sets the parameter."
                                            if(length(value) > 1 || value == 0){stop("You can't enter a zero value.")}
                                            SDstringency <<- value
                                          },
                                        
                                        setSettings =
                                          function(...){
                                            settings <- list(...)
                                            parameters <- names(settings)
                                            for(i in 1:length(settings)){
                                              if(parameters[i] == "ManualThresholds"){
                                                setManualThresholds(settings[[i]])
                                              }
                                              if(parameters[i] == "AutoThresholds"){
                                                AutoThresholds <<- settings[[i]]
                                              }
                                              if(parameters[i] == "ManualFallback"){
                                                ManualFallback <<- settings[[i]]
                                              }
                                              if(parameters[i] == "SDstringency"){
                                                setSDstringency(settings[[i]])
                                              }
                                            }
                                          },
                                        
                                        textSummary =
                                          function(){
                                            return(paste0('Settings for detecting blocks from recombination signal:\n',
                                                          '--------------------------------------------------------\n',
                                                          'Manual sequence similarity thresholds (ManualThresholds): ',
                                                          paste(ManualThresholds, collapse=", "),
                                                          '\n\nAutomatically decide on thresholds (AutoThresholds): ',
                                                          AutoThresholds,
                                                          '\n\nFall back to manual thresholds (ManualFallback): ',
                                                          ManualFallback,
                                                          '\n\nStandard deviation divisor (SDStringency): ', SDstringency
                                            ))
                                          },
                                        
                                        show =
                                          function(){
                                            cat(textSummary())
                                          }
                                      )
)

#' A reference class storing the settings for recombination block dating. 
#' @name BlockDatingSettings
#' @field MutationRate Numeric vector of length one. Stores the mutation rate to be used when dating blocks.
#' @field PValue Numeric vector of length one. Stores the critical alpha value for testing the signifcance of recombination regions.
#' @field BonfCorrection Logical vector of length one, stores the option of whether the critical value stored in PValue will be corrected.
#' @field DateAnyway Logical vector of length one, sotres the option of whether blocks will be dated despite failing the critical alpha.
#' @field MutationCorrection Character vector of length 1, can be any of the model supported by the ape package. Default is "JC69". 
BlockDatingSettings <- setRefClass("BlockDatingSettings",
                                   
                                   fields = list(
                                     MutationRate = "numeric",
                                     PValue = "numeric",
                                     BonfCorrection = "logical",
                                     DateAnyway = "logical",
                                     MutationCorrection = "character"
                                   ),
                                   
                                   methods = list(
                                     initialize =
                                       function(){
                                         MutationRate <<- 10e-09
                                         PValue <<- 0.005
                                         BonfCorrection <<- TRUE
                                         DateAnyway <<- FALSE
                                         MutationCorrection <<- "JC69"
                                       },
                                     
                                     setMutationRate =
                                       function(newRate){
                                         "Sets a new mutation rate for the settings."
                                         if(length(newRate) > 1){stop("Input must be a single value.")}
                                         MutationRate <<- newRate
                                       },
                                     
                                     setPValue =
                                       function(newValue){
                                         "Set a new critical value for significance testing of recombinant blocks."
                                         if(length(newValue) > 1){stop("Input must be a single value.")}
                                         PValue <<- newValue
                                       },
                                     
                                     setBonferonni =
                                       function(newBonf){
                                         "Set whether the critical value should be subject to bonferroni correction during block significance testing."
                                         if(length(newBonf) > 1){stop("Input must be a single value.")}
                                         BonfCorrection <<- newBonf
                                       },
                                     
                                     setDateAnyway =
                                       function(newValue){
                                         "Set whether blocks should be kept and dated even if they fail the significance test."
                                         if(length(newValue) > 1){stop("Input must be a single value.")}
                                         DateAnyway <<- newValue
                                       },
                                     
                                     setMutationCorrection =
                                       function(model){
                                         "Set the model of sequence evolution to correct the distances/number of mutations used in block dating algorithm."
                                         if(length(model) > 1){stop("Input must be a single value.")}
                                         if(!any(model == c("raw", "TS", "TV", "JC69", "K80", "F81",
                                                            "K81", "F84", "BH87", "T92", "TN93", "GG95"))){
                                           stop(paste0("Provided model must be one of the following: ", paste(c("raw", "TS", "TV", "JC69", "K80", "F81",
                                                                                                                "K81", "F84", "BH87", "T92", "TN93", "GG95."), collapse=", ")))
                                         }
                                         MutationCorrection <<- model
                                       },
                                     
                                     setSettings =
                                       function(...){
                                         settings <- list(...)
                                         parameters <- names(settings)
                                         for(i in 1:length(settings)){
                                           if(parameters[i] == "MutationCorrection"){
                                             setMutationCorrection(settings[[i]])
                                           }
                                           if(parameters[i] == "DateAnyway"){
                                             setDateAnyway(settings[[i]])
                                           }
                                           if(parameters[i] == "BonfCorrection"){
                                             setBonferonni(settings[[i]])
                                           }
                                           if(parameters[i] == "PValue"){
                                             setPValue(settings[[i]])
                                           }
                                           if(parameters[i] == "MutationRate"){
                                             setMutationRate(settings[[i]])
                                           }
                                         }
                                       },
                                     
                                     textSummary =
                                       function(){
                                         return(paste0('Settings for testing and dating recombination blocks:\n',
                                                       '-----------------------------------------------------\n',
                                                       'Assumed substitution rate for dating (MutationRate): ',
                                                       MutationRate,
                                                       '\n\nCritical alpha value for significance testing (PValue): ',
                                                       PValue,
                                                       '\n\nApply bonferroni correction to critical alpha (BonfCorrection): ',
                                                       BonfCorrection,
                                                       '\n\nKeep and date blocks that fail the alpha (DateAnyway): ', DateAnyway,
                                                       '\n\nAssumed mutation model for dating (MutationCorrection): ', MutationCorrection))
                                       },
                                     
                                     show =
                                       function(){
                                         cat(textSummary())
                                       }
                                   ))