library(shinydashboard)
library(HybridCheck)

options(shiny.maxRequestSize=10000*1024^2)

comboChoiceSorter <- function(inputs){
  return(lapply(inputs, function(x){
    if(x == "NOT.TESTED" || x == "TESTED" || x == "ALL"){
      return(x)
    } else {
      chars <- unlist(
        strsplit(
          unlist(strsplit(x, split = ", ")),
          ": "))
      return(c(chars[2], chars[4], chars[6], chars[8]))
    }
  }))
}

function(input, output, session){
  
  HCobj <- HC$new()
  
  updateSequence <- reactive({
    if(!is.null(input$fastafile$datapath) && (grepl(".fas", input$fastafile$name) || grepl(".fasta", input$fastafile$name)))
      HCobj$inputDNA(input$fastafile$datapath)
  })
  
  updatePopulations <- reactive({
    input$setPops
    if(HCobj$DNA$hasDNA()){
      if(!isolate(input$oneSeqOnePop)){
        HCobj$setPopulations(NULL)
      } else {
        nameSel <- unlist(lapply(1:isolate(input$numPops), function(i) paste0("PopulationName", i)))
        popSel <- unlist(lapply(1:isolate(input$numPops), function(i) paste0("Population", i)))
        populations <- setNames(object = lapply(popSel, function(i) isolate(input[[i]])), 
                                unlist(lapply(nameSel, function(i) isolate(input[[i]]))))
        HCobj$setPopulations(populations)
      }
    }
  })
  
  output$seqNum <- renderText({
    updateSequence()
    validate(need(HCobj$DNA$hasDNA(), "0"))
    HCobj$DNA$numberOfSequences()
  })
  
  output$bpLen <- renderText({
    updateSequence()
    validate(need(HCobj$DNA$hasDNA(), "0"))
    HCobj$DNA$getFullLength()
  })
  
  output$infLen <- renderText({
    updateSequence()
    validate(need(HCobj$DNA$hasDNA(), "0"))
    HCobj$DNA$getInformativeLength()
  })
  
  output$numPopulations <- renderText({
    updateSequence()
    updatePopulations()
    validate(need(HCobj$DNA$hasDNA(), "0"))
    validate(need(HCobj$DNA$hasPopulations(), "0"))
    HCobj$DNA$numberOfPopulations()
  })
  
  output$seqNames <- renderText({
    updateSequence()
    validate(need(HCobj$DNA$hasDNA(), "Sequences have not been loaded."))
    paste(HCobj$DNA$getSequenceNames(), collapse = ", ")
  })
  
  output$popDetails <- renderUI({
    updateSequence()
    updatePopulations()
    popNames <- HCobj$DNA$namesOfPopulations()
    validate(need(HCobj$DNA$hasPopulations(), "No populations defined."))
    lapply(1:HCobj$DNA$numberOfPopulations(), function(i){
      fluidRow(box(title = popNames[i], solidHeader = TRUE, width = 12,
                   paste0(HCobj$DNA$Populations[[i]], collapse = ", ")
      ))
    }) 
  })
  
  output$populationsGen <- renderUI({
    updateSequence()
    validate(need(HCobj$DNA$hasDNA(), "Sequences have not been loaded."))
    validate(need(!is.na(input$numPops), "Enter a number of populations."))
    lapply(1:input$numPops, function(i) {
      fluidRow(
        column(6, textInput(inputId=paste0("PopulationName", i), label=paste0("Population Name"))),
        column(6,
               selectInput(inputId=paste0("Population", i), label=paste0("Sequences"),
                           choices = HCobj$DNA$getSequenceNames(), multiple = TRUE) 
        )
      )
    })
  })
  
  output$ABBAtree <- renderImage({
    width  <- session$clientData$output_ABBAtree_width
    height <- session$clientData$output_ABBAtree_height
    pixelratio <- session$clientData$pixelratio
    outfile <- tempfile(fileext='.png')
    png(outfile, width = width * pixelratio, height = height * pixelratio,
        res = 72 * pixelratio)
    plot(ape:::read.tree(text="(((P1,P2),P3),P4);"))
    ape:::nodelabels("A", c(1,4), adj = c(-2.5, 0.5), bg = "red", col="white")
    ape:::nodelabels("B", c(2,3), adj = c(-2.5, 0.5), bg = "blue", col="white")
    dev.off()
    list(src = outfile,
         width = width,
         height = height,
         alt = "This is alternate text")}, deleteFile = TRUE)
  
  output$BABAtree <- renderImage({
    width  <- session$clientData$output_BABAtree_width
    height <- session$clientData$output_BABAtree_height
    pixelratio <- session$clientData$pixelratio
    outfile <- tempfile(fileext='.png')
    png(outfile, width=width*pixelratio, height=height*pixelratio,
        res=72*pixelratio)
    plot(ape:::read.tree(text="(((P1,P2),P3),P4);"))
    ape:::nodelabels("A", c(2,4), adj = c(-2.5, 0.5), bg = "red", col="white")
    ape:::nodelabels("B", c(1,3), adj = c(-2.5, 0.5), bg = "blue", col="white")
    dev.off()
    list(src = outfile,
         width = width,
         height = height,
         alt = "This is alternate text")}, deleteFile = TRUE)
  
  output$fttGen <- renderUI({
    updateSequence()
    updatePopulations()
    validate(need(!is.na(input$fttNumCombos), "Enter a number of taxa combos to analyse."))
    validate(need(length(HCobj$DNA$Populations) >= 4, "You need 4 or more populations defined."))
    fluidRow(
      column(3,
             lapply(1:input$fttNumCombos, function(i) {
               selectInput(inputId=paste0("P1", i), label = "P1",
                           choices = HCobj$DNA$namesOfPopulations())
             })),
      column(3,
             lapply(1:input$fttNumCombos, function(i) {
               selectInput(inputId=paste0("P2", i), label = "P2",
                           choices = HCobj$DNA$namesOfPopulations())
             })),
      column(3,
             lapply(1:input$fttNumCombos, function(i) {
               selectInput(inputId=paste0("P3", i), label = "P3",
                           choices = HCobj$DNA$namesOfPopulations())
             })),
      column(3,
             lapply(1:input$fttNumCombos, function(i) {
               selectInput(inputId=paste0("A", i), label = "P4",
                           choices = HCobj$DNA$namesOfPopulations())
             }))
    )
  })
  
  updateTaxaCombos <- reactive({
    input$setCombinations
    if(isolate(input$fttAutoSets)){
      HCobj$prepareFourTaxonTests()
    } else {
      combos <- lapply(1:isolate(input$fttNumCombos), function(i){
        c(P1 = isolate(input[[paste0("P1", i)]]), P2 = isolate(input[[paste0("P2", i)]]),
          P3 = isolate(input[[paste0("P3", i)]]), A = isolate(input[[paste0("A", i)]]))
      })
      combos <- combos[unlist(lapply(combos, function(x) length(unique(x)) == 4))]
      if(length(combos) > 0)
        HCobj$prepareFourTaxonTests(combos)
    }
  })
  
  output$generatedCombos <- renderUI({
    updateTaxaCombos()
    lapply(HCobj$FTTmodule$results, function(x){
      box(title = NULL, width = 12,
          paste0("P1: ", x$P1, ",\nP2: ", x$P2, ",\nP3: ", x$P3, ",\nA: ", x$A))
    })
  })
  
  output$combosToAnalyze <- renderUI({
    updateTaxaCombos()
    validate(need(HCobj$FTTmodule$hasTaxaCombos(), "Set taxa combinations to test."))
    column(width = 12,
           selectInput("fttToAnalyze", 
                       tags$strong("Run / rerun four taxon test for combination:"),
                       c("NOT.TESTED", "ALL", "TESTED", HCobj$FTTmodule$printAllNames()),
                       selected = "ALL", multiple = TRUE),
           numericInput("fttBlockSize", "Block Size:", NULL),
           numericInput("fttNumberOfBlocks", "Number of Blocks:", NULL),
           actionButton("runFTTs", "Run Tests")
    )
  })
  
  doFTTests <- reactive({
    input$runFTTs
    if((length(HCobj$FTTmodule$results) != 0) && (input$runFTTs > 0)){
      selections <- comboChoiceSorter(isolate(input$fttToAnalyze))
      HCobj$runFourTaxonTests(selections = selections,
                                   numberOfBlocks = isolate(input$fttNumberOfBlocks),
                                   blockLength = isolate(input$fttBlockSize))
    }  
  })
  
  output$combosToView <- renderUI({
    updateTaxaCombos()
    validate(need(HCobj$FTTmodule$hasTaxaCombos(), "Set taxa combinations to test."))
    column(width = 12,
           selectInput("fttView", 
                       tags$strong("View four taxon test results for combination:"),
                       c("NOT.TESTED", "ALL", "TESTED", "SIGNIFICANT", HCobj$FTTmodule$printAllNames()),
                       selected = "ALL", multiple = TRUE)
    )
  })
  
  output$fttResults <- renderUI({
    doFTTests()
    selectedTests <- HCobj$FTTmodule$getFTTs(comboChoiceSorter(input$fttView))
    lapply(selectedTests, function(x){
      if(!x$noTestPerformed()){
        testName <- paste0("P1: ", x$P1, ", P2: ", x$P2, ", P3: ", x$P3,
                           ", P4: ", x$A)
        tableName <- paste0("table_", x$P1, x$P2, x$P3, x$A)
        table <- x$table[, c("BlockStart", "BlockEnd",
                              "S", "numBinomialP",
                              "ABBA", "BABA", "D", "Fd_1DD4_D0", "Fd_D2D4_D0",
                              "blockFraction", "scaledPseudoD", "scaledPseudoFd_1DD4",
                              "scaledPseudoFd_D2D4")]
        colnames(table) <- c("BlockStart", "BlockEnd", "S", "BinomP", "ABBA", "BABA",
                             "D", "Fd:P2,P3", "Fd:P1,P3", "BlockFraction", "PseudoD",
                             "PseudoFd:P2,P3", "PeudoFd:P1,P3")
        table[is.na(table)] <- 0
        output[[tableName]] <- renderTable(table)
        box(title = testName, width = 12, solidHeader = TRUE,
            status = "primary", 
            fluidRow(
              box(title = "Global Stats", status = "primary", width = 12,
                  fluidRow(
                    valueBox(round(x$ABBA, 2), "ABBA", width = 3, color = "orange"),
                    valueBox(round(x$BABA, 2), "BABA", width = 3, color = "yellow"),
                    valueBox(round(x$X2_P, 10), "P-Value", width = 3, color = "green"),
                    valueBox(round(x$Observed_D, 5), "D", width = 3, color = "blue"),
                    valueBox(round(x$D_jEstimate, 5), "D (Jackknifed)", width = 3, color = "blue"),
                    valueBox(round(x$Observed_Fd_1DD4, 5), "Fd(P2, P3)", width = 3, color = "navy"),
                    valueBox(round(x$Fd_1DD4_jEstimate, 5), "Fd(P2, P3) (Jackknifed)", width = 3, color = "navy"),
                    valueBox(round(x$Observed_Fd_D2D4, 5), "Fd(P1, P3)", width = 3, color = "navy"),
                    valueBox(round(x$Fd_D2D4_jEstimate, 5), "Fd(P1, P3) (Jackknifed)", width = 3, color = "navy"),
                    valueBox(round(x$D_jZ, 5), "Z score for D", width = 3, color = "purple"),
                    valueBox(round(x$Fd_1DD4_jZ, 5), "Z score for Fd(P2, P3)", width = 3, color = "purple"),
                    valueBox(round(x$Fd_D2D4_jZ, 5), "Z score for Fd(P1, P3)", width = 3, color = "purple")
                  )
              ),
              box(title = "Jack-knifed segment stats", status = "warning", width = 12,
                  tableOutput(tableName)
              )
            )
        )
      }
    })
  })
  
  updateTripletGenSettings <- reactive({
    input$comboGen
    newoptions <- c()
    if(isolate(input$fttOrPops) == "ftt"){
      newoptions <- c(newoptions, 1L)
    } else {
      newoptions <- c(newoptions, 2L)
    }
    if(isolate(input$distancebased)){
      if(isolate(input$distancemethod) == "yes"){
        newoptions <- c(newoptions, 4L)
      } else {
        newoptions <- c(newoptions, 3L)
      }
    }
    validate(
      need(!is.null(newoptions), "Select either to generate triplet combinations based on group specification or based on distance, or both.")
    )
    validate(
      need(HCobj$DNA$hasDNA(), "No sequence file is loaded into HC.")
    )
    HCobj$setParameters("TripletGeneration",
                             Method = newoptions,
                             DistanceThreshold = isolate(input$mandistthreshold),
                             PartitionStrictness = as.integer(isolate(input$partitionStrictness)))
  })
  
  output$NumCombos <- renderText({
    updateTripletGenSettings()
    length(HCobj$comparrisonSettings$AcceptedCombinations)
  })
  
  output$GeneratedTriplets <- renderText({
    updateTripletGenSettings()
    HCobj$comparrisonSettings$htmlCombinations()
  })
  
  output$TripletToAnalyze <- renderUI({
    updateTripletGenSettings()
    selectInput("tripletToAnalyze", 
                tags$strong("Run / rerun analysis for triplets:"),
                c("ALL", HCobj$comparrisonSettings$printAcceptedCombinations()),
                selected = "ALL", multiple = TRUE)
  })
  
  updatePlottingSettings <- reactive({
    HCobj$setParameters("Plotting", PlotTitle = input$plotTitle,
                             TitleSize = input$plotTitleSize,
                             TitleFace = input$plotTitleFace,
                             TitleColour = input$plotTitleColour,
                             XLabels = input$plotXLabels,
                             YLabels = input$plotYLabels,
                             XTitle = input$plotXTitle,
                             XTitleFontSize = input$plotXTitleFontSize,
                             XTitleColour = input$plotXTitleColour,
                             XLabelSize = input$plotXLabelSize,
                             XLabelColour = input$plotXLabelColour,
                             YTitle = input$plotYTitle,
                             YTitleFontSize = input$plotYTitleFontSize,
                             YTitleColour = input$plotYTitleColour,
                             YLabelSize = input$plotYLabelSize,
                             YLabelColour = input$plotYLabelColour,
                             Legends = input$plotLegends,
                             LegendFontSize = input$plotLegendFontSize,
                             MosaicScale = input$plotMosaicScale)
  })
  
  analysis <- reactive({
    input$analysisGO
    if(HCobj$DNA$hasDNA() && (input$analysisGO != 0)){
      HCobj$setParameters("SSAnalysis", WindowSize = as.integer(isolate(input$windowSize)),
                               StepSize = as.integer(isolate(input$stepSize)))
      HCobj$setParameters("BlockDetection", ManualThresholds = isolate(input$manBlockDetectDist),
                               AutoThresholds = (isolate(input$detectionMethod) == "no"),
                               ManualFallback = isolate(input$fallbackManual))
      dateanyway <- !isolate(input$eliminateinsignificant)
      HCobj$setParameters("BlockDating", MutationRate = isolate(input$mu), PValue = isolate(input$alpha),
                               BonfCorrection = isolate(input$bonf), DateAnyway = dateanyway,
                               MutationCorrection = isolate(input$correctionModel))
      selections <- strsplit(isolate(input$tripletToAnalyze), ", ")
      HCobj$triplets$generateTriplets(HCobj$DNA, HCobj$comparrisonSettings, HCobj$filesDirectory)
      tripletsToDo <- HCobj$triplets$getTriplets(selections)
      numToDo <- length(tripletsToDo)
      prog <- Progress$new(session, min = 1, max = numToDo)
      prog$set(value = 0, message = "Scanning SS in triplets...")
      for(i in 1:numToDo){
        HybridCheck:::scan.similarity(HCobj$DNA, tripletsToDo[[i]], HCobj$ssAnalysisSettings)
        prog$set(value = i)
      }
      prog$close()
      prog <- Progress$new(session, min = 1, max = numToDo)
      prog$set(value = 0, message = "Finding blocks in triplets...")
      for(i in 1:numToDo){
        tripletsToDo[[i]]$putativeBlockFind(HCobj$blockDetectionSettings)
        prog$set(value = i)
      }
      prog$close()
      prog <- Progress$new(session, min = 1, max = numToDo)
      prog$set(value = 0, message = "Dating blocks...")
      for(i in 1:numToDo){
        tripletsToDo[[i]]$blockDate(HCobj$DNA, HCobj$blockDatingSettings)
        prog$set(value = i)
      }
      prog$close() 
    }
  })
  
  output$TripletSelector <- renderUI({
    analysis()
    selectInput("tripletSelection", tags$strong("View triplet:"),
                HCobj$triplets$tripletCombinations())
  })
  
  output$barsPlot <- renderPlot({
    analysis()
    updatePlottingSettings()
    validate(need(is.character(input$tripletSelection), 
                  "Triplets to choose from have not been generated yet."))
    selection <- unlist(strsplit(input$tripletSelection, ", "))
    validate(need(!HCobj$triplets$getTriplets(selection)[[1]]$noScanPerformed(),
                  "Sequence similarity scan has not been performed for this triplet yet."))
    barsPlot <- HCobj$plotTriplets(unlist(strsplit(input$tripletSelection, ", ")), What="Bars")[[1]]
    print(barsPlot)
  })
  
  output$linesPlot <- renderPlot({
    analysis()
    updatePlottingSettings()
    validate(need(is.character(input$tripletSelection), 
                  "Triplets to choose from have not been generated yet."))
    selection <- unlist(strsplit(input$tripletSelection, ", "))
    validate(need(!HCobj$triplets$getTriplets(selection)[[1]]$noScanPerformed(),
                  "Sequence similarity scan has not been performed for this triplet yet."))
    linesPlot <- HCobj$plotTriplets(selection, What="Lines")[[1]]
    print(linesPlot)
  })
  
  output$blocksTable <- renderDataTable({
    analysis()
    validate(need(is.character(input$tripletSelection),
                  "Triplets to choose from have not been generated yet."))
    validate(
      need(!HCobj$triplets$getTriplets(
        unlist(strsplit(input$tripletSelection, ", ")))[[1]]$blocksNotDated(),
        "Recombinant regions have not been found or dated yet for the selected triplet."))
    HCobj$tabulateDetectedBlocks(unlist(
      strsplit(input$tripletSelection, ", ")), Neat=TRUE)[, c(-1, -3, -6, -7)]
  })
  
  output$saveFTT <- downloadHandler(filename = 
                                        function(){
                                          paste0(strsplit(input$fastafile$name, ".fas")[1], "_FTT.csv")
                                        },
                                      content =
                                        function(file){
                                          write.csv(HCobj$tabulateFourTaxonTests(unlist(strsplit(input$combosToView, ", "))), file)
                                        }
  )
  
  output$saveTable <- downloadHandler(filename = 
                                        function(){
                                          paste0(strsplit(input$fastafile$name, ".fas")[1], "_Triplet_",
                                                 paste(unlist(strsplit(input$tripletSelection, ", ")), collapse = ":"), ".csv")
                                        },
                                      content =
                                        function(file){
                                          write.csv(HCobj$tabulateDetectedBlocks(unlist(strsplit(input$tripletSelection, ", ")), Neat=TRUE), file)
                                        }
  )
  
  output$saveBarPlots <- downloadHandler(filename =
                                           function(){
                                             paste0(strsplit(input$fastafile$name, ".fas")[1], "_Triplet_",
                                                    paste(unlist(strsplit(input$tripletSelection, ", ")), collapse = ":"), "_Bars.png")
                                           },
                                         content = 
                                           function(file){
                                             selection <- unlist(strsplit(isolate(input$tripletSelection), ", "))
                                             whichToPlot <- isolate(input$whichPlotToSave)
                                             if(whichToPlot == "Both"){
                                               whichToPlot <- c("Bars", "Lines")
                                             }
                                             Plot <- HCobj$plotTriplets(selection, What = whichToPlot)[[1]]
                                             ggsave("plot.png", plot = Plot, width = isolate(input$widthSave), height = isolate(input$heightSave),
                                                    dpi = isolate(input$resSave), units = "in")
                                             file.copy("plot.png", file, overwrite=TRUE)
                                           })
  
  output$userBlocksSeqSelect <- renderUI({
    updateSequence()
    selectInput("seqChoice", "Select two sequences",
                HCobj$DNA$getSequenceNames(), multiple = TRUE)
  })
  
  addUserBlocks <- reactive({
    input$addUBButton
    validate(need(isolate(input$startPosition) < isolate(input$endPosition), "Start position must be less than the end position."))
    if(length(isolate(input$seqChoice)) == 2){
      HCobj$addUserBlock(isolate(input$seqChoice), isolate(input$startPosition), isolate(input$endPosition))
    }
  })
  
  clearUserBlocks <- reactive({
    input$clearUBButton
    if(length(isolate(input$seqChoice)) == 2){
      HCobj$clearUserBlocks(isolate(input$seqChoice))
    }
  })
  
  dateUserBlocks <- reactive({
    input$dateUBButton
    HCobj$setParameters("BlockDating", MutationRate = isolate(input$mu2), PValue = isolate(input$alpha2),
                             BonfCorrection = isolate(input$bonf2), DateAnyway = !isolate(input$eliminateinsignificant2),
                             MutationCorrection = isolate(input$correctionModel2))
    HCobj$dateUserBlocks()
  })
  
  output$userBlocksTable <- renderDataTable({
    updateSequence()
    clearUserBlocks()
    addUserBlocks()
    dateUserBlocks()
    HCobj$tabulateUserBlocks()
  })
  
}
