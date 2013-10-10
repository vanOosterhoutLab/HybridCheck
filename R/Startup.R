.onAttach <- function(...) {
  cat("\nTACTTTGTACCTAAGTATGCATTACGTTACGTTAGTAGCTGGACCTAGTAAATCGGA     
,--.  ,--.         ,--.   ,------. ,--.,------.   ,---.
|  '--'  |,--. ,--.|  |-. |  .--. '|  ||  .-.  \\ '   .-'
|  .--.  | \\  '  / | .-. '|  '--'.'|  ||  |  \\  :`.  `-.
|  |  |  |  \\   '  | `-' ||  |\\  \\ |  ||  '--'  /.-'    |
`--'  `--'.-'  /    `---' `--' '--'`--'`-------' `-----'
          `---'
ATGAAACATGGATTCATACG - Version 1.0 - CGACCTGGATCATTTAGCCT\n\n
Hybridisation, Recombination and Introgression Detection
                   and Dating Package.\n
-----------------=========****=========------------------
Cite: TBD
Licence: GPL (Like R and most packages).
http://ward9250.github.io/HybRIDS/
-----------------=========****=========------------------\n")
  cat("\nEnter HybRIDS_GUI() at the console prompt for the interactive GUI.\n\n")
}

#' Start a GUI for HybRIDS
#' 
#' @export
HybRIDS_GUI <- function(toolkit="tcltk"){
  require(gWidgets)
  require(gWidgetstcltk)
  options(guiToolkit=toolkit)
  invisible(HybRIDS_gui_generator$new())
}

HybRIDS_gui_generator <- setRefClass(
  "HybRIDSGui",
  fields = list(
    #Widgets
    MainWindow = "ANY",
    DataButton = "ANY",
    TripletsButton = "ANY",
    SSButton = "ANY",
    BlockDetectButton = "ANY",
    BlockDatingButton = "ANY",
    PlottingButton = "ANY",
    AnalysisSwitcher = "ANY",
    SecondMainwinGroup = "ANY",
    DataSelection = "ANY",
    DataSelectGroup2 = "ANY",
    SwitcherButton = "ANY",
    DataSelectGroup = "ANY",
    TripletFilter = "ANY",
    DataSummary = "ANY",
    DataPlotFrame = "ANY",
    AnalysisFrame = "ANY",
    StatusBar = "ANY",
    # Data
    HybRIDS_Sessions = "list",
    CommandHistory = "list",
    CurrentPlot = "ANY"
  ),
  methods = list(
    initialize = function()
    {
      ## Initializing all of the GUI components...
      MainWindow <<- gwindow("HybRIDS Control Panel", visible = FALSE)
      DataSelection <<- gframe("Data Selector", cont=MainWindow, horizontal=F)
      DataSelectGroup2 <<- gframe("HybRIDS Session Selection", cont=DataSelection, horizontal=T)
      AnalysisSwitcher <<- gcombobox("", cont=DataSelectGroup2, width=100, editable=T)
      SwitcherButton <<- gbutton("Delete Analysis", cont=DataSelectGroup2, handler=function(h,...){DelSession(svalue(AnalysisSwitcher))})
      DataSelectGroup <<- gframe("Triplet Selection Filter", cont=DataSelection, horizontal=T)
      TripletFilter <<- gedit("ALL", width=20, cont=DataSelectGroup)
      AnalysisFrame <<- gframe("Analysis", cont=MainWindow, horizontal=F, handler=function(h,...){SSHandler()})
      DataButton <<- gbutton("Step 1). Select and Load DNA Data (FASTA)", container=AnalysisFrame, expand=TRUE, handler=function(h,...){LoadIn()})
      TripletsButton <<- gbutton("Step 2). Generate the Triplet Sets for SSAnalysis", container=AnalysisFrame, expand=TRUE, handler=function(h,...){TripletsHandler()})
      SSButton <<- gbutton("Step 3). Sequence Similarity Analysis", container=AnalysisFrame, expand=TRUE, handler=function(h,...){SSHandler()})
      BlockDetectButton <<- gbutton("Step 4). Detect Putative Blocks", container=AnalysisFrame, expand=TRUE, handler=function(h,...){BlockDetectHandler()})
      BlockDatingButton <<- gbutton("Step 5). Test and Date Putative Blocks", container=AnalysisFrame, expand=TRUE, handler=function(h,...){TestDateHandler()})
      DataPlotFrame <<- gframe("Data and Plotting", cont=MainWindow)
      PlottingButton <<- gbutton("Data Export & Plotting", container=DataPlotFrame, expand=TRUE, handler=function(h,...){PlotExportHandler()})
      DataSummary <<- gbutton("HybRIDS Session Summary", cont=DataPlotFrame, expand=T,
                              handler = function(h, ...){
                                AnalysisDisplay(svalue(AnalysisSwitcher))
                              })
      StatusBar <<- gstatusbar("Ready...", container=MainWindow)
      visible(MainWindow) <<- TRUE
    },
    TripletsHandler =
      function(){
        sel <- svalue(AnalysisSwitcher)
        TripletWindow <- gwindow(paste("GenerateTriplets for Session:", sel))
        preanalysis1 <- gcheckbox("Pre-Sort by Raw Sequence Similarity Threshold", checked = HybRIDS_Sessions[[sel]]$TripletParams$Method==2, use.togglebutton=FALSE, handler = function(h,...){
          if(svalue(preanalysis1)==TRUE){
            if(svalue(preanalysis2)==TRUE){
              svalue(preanalysis2) <- FALSE
            }
            HybRIDS_Sessions[[sel]]$TripletParams$Method <<- 2
          } else {
            HybRIDS_Sessions[[sel]]$TripletParams$Method <<- 1
          }
        }, action = NULL, container = TripletWindow)
        labpre <- glabel("Pre-Sort Raw Sequence Similarity Threshold", cont=TripletWindow)
        prethreshold <- gedit(text=as.character(HybRIDS_Sessions[[sel]]$TripletParams$SortThreshold), width = 20, cont=TripletWindow)
        preanalysis2 <- gcheckbox("Pre-Sort by Distance Distribution Dependent Cutoff", checked = HybRIDS_Sessions[[sel]]$TripletParams$Method==3, use.togglebutton=FALSE, handler = function(h,...){
          if(svalue(preanalysis2)==TRUE){
            if(svalue(preanalysis1)==TRUE){
              svalue(preanalysis1) <- FALSE
            }
            HybRIDS_Sessions[[sel]]$TripletParams$Method <<- 3
          } else {
            HybRIDS_Sessions[[sel]]$TripletParams$Method <<- 1
          }
        }, action = NULL, container = TripletWindow)
        tripletsetbutton <- gbutton("Generate Triplets", cont=TripletWindow, handler=function(h,...){
          HybRIDS_Sessions[[sel]]$TripletParams$SortThreshold <<- as.numeric(svalue(prethreshold))
          svalue(StatusBar) <<- paste("Generating Triplets for", sel)
          HybRIDS_Sessions[[sel]]$makeTripletCombos()
          gmessage("Made all the triplets", icon="info")
          svalue(StatusBar) <<- "Ready..."
          dispose(TripletWindow)
        })
      },
    AnalysisDisplay =
      function(sel) {
        CheckPresent(sel)
        t <- capture.output(HybRIDS_Sessions[[sel]]$show())
        tempwin <- gwindow(title=paste(sel,":", t[1], sep=" "), cont=F, visible=F)
        f1 <- gframe(paste(t[3]), cont=tempwin)
        l1 <- glabel(paste(t[7:9], sep="\n"), cont=f1)
        f2 <- gframe(paste(t[11]), cont=tempwin)
        l2 <- glabel(paste(t[c(13,14)], sep="\n"), cont=f2)
        f3 <- gframe(paste(t[16]), cont=tempwin)
        l3 <- glabel(paste(t[18:19], sep="\n"), cont=f3)
        f4 <- gframe(paste(t[21]), cont=tempwin)
        l4 <- glabel(paste(t[23:25], sep="\n"), cont=f4)
        f5 <- gframe(paste(t[27]), cont=tempwin)
        l5 <- glabel(paste(t[29:30], sep="\n"), cont=f5)
        l6 <- glabel(t[32], cont=tempwin)
        but1 <- gbutton("Print as text to console", cont=tempwin, handler = function(h,...){HybRIDS_Sessions[[sel]]$show()})
        visible(tempwin) <- T
      },
    LoadIn =
      function(){
        sequencesname <- ginput(message="Give a name to the HybRIDS session you're creating", title="Sequence Names")
        filename <- gfile(text="Choose a file", 
                          type="open",
                          filter = list("Fasta Files" = list(patterns = c("*.fasta","*.fas"))),
                          action="print",
                          handler= 
                            function(h,...){
                              do.call(h$action, list(h$file))
                            }
        )
        svalue(StatusBar) <<- "Loading Data & Creating HybRIDS session..."
        HybRIDS_Sessions <<- append(HybRIDS_Sessions, HybRIDS$new(filename))
        names(HybRIDS_Sessions)[[length(HybRIDS_Sessions)]] <<- sequencesname
        AnalysisSwitcher[] <<- names(HybRIDS_Sessions)
        gmessage("Loaded DNA Data and created new HybRIDS session", icon="info")
        svalue(StatusBar) <<- "Ready..."
      },
    DelSession =
      function(sel){
        HybRIDS_Sessions[which(names(HybRIDS_Sessions) %in% sel)] <<- NULL
        AnalysisSwitcher[] <<- names(HybRIDS_Sessions)
        if(length(HybRIDS_Sessions) > 0){
          svalue(AnalysisSwitcher) <<- names(HybRIDS_Sessions[[1]])
        } else {
          svalue(AnalysisSwitcher) <<- ""
        }
      },
    SSHandler =
      function(){
        sel <- svalue(AnalysisSwitcher)
        SSWindow <- gwindow( paste("Doing SSAnalysis for", sel), visible=F)
        subwinopt <- gframe("Window Settings", cont = SSWindow, horizontal=F)
        lab3 <- glabel("Window Size",cont=subwinopt)
        winsizespin <- gedit(text=as.character(HybRIDS_Sessions[[svalue(AnalysisSwitcher)]]$SSAnalysisParams$WindowSize), width = 20, cont=subwinopt)
        lab4 <- glabel("Step Size",cont=subwinopt)
        stepsizespin <- gedit(text=as.character(HybRIDS_Sessions[[svalue(AnalysisSwitcher)]]$SSAnalysisParams$StepSize), width = 20, cont=subwinopt)
        TripletFilterCheck <- gcheckbox("Apply Triplet Selection Filter",
                                        cont=SSWindow)
        Analyze1 <- gbutton("Analyze", container=SSWindow, handler=function(h,...){
          HybRIDS_Sessions[[sel]]$SSAnalysisParams$WindowSize <<- as.numeric(svalue(winsizespin))
          HybRIDS_Sessions[[sel]]$SSAnalysisParams$StepSize <<- as.numeric(svalue(stepsizespin))
          if(svalue(TripletFilterCheck)==TRUE){
            HybRIDS_Sessions[[sel]]$analyzeSS(svalue(TripletFilter))
          } else {
            HybRIDS_Sessions[[sel]]$analyzeSS()
          }
          gmessage("Finished analyzing sequence similarities of Sequence Triplets", icon="info")
          dispose(SSWindow)
        })
        visible(SSWindow) <- TRUE
      },
    BlockDetectHandler =
      function(){
        sel <- svalue(AnalysisSwitcher)
        BlockDetectWindow <- gwindow(paste("Block Detection for HybRIDS session:", sel), visible=FALSE)
        ThresholdsLabel <- glabel("Manual Similarity Thresholds (%)", cont=BlockDetectWindow)
        ManualThresholdsEntry <- gedit(paste(as.character(HybRIDS_Sessions[[sel]]$BlockDetectionParams$ManualThresholds), collapse=","), width=10, cont=BlockDetectWindow)
        AutoDetectCheck <- gcheckbox("Auto-detect Thresholds?", checked = HybRIDS_Sessions[[sel]]$BlockDetectionParams$AutoThresholds, use.togglebutton = FALSE, cont=BlockDetectWindow)
        ManualFallbackCheck <- gcheckbox("Fallback to manual values?", checked = HybRIDS_Sessions[[sel]]$BlockDetectionParams$ManualFallback, use.togglebutton = FALSE, cont=BlockDetectWindow)
        FilterCheck <- gcheckbox("Apply Triplet Filter", cont=BlockDetectWindow)
        SDstringencyLab <- glabel("SD Stringency", cont=BlockDetectWindow)
        SDstringencyGedit <- gedit(as.character(HybRIDS_Sessions[[sel]]$BlockDetectionParams$SDstringency), width=5, cont=BlockDetectWindow)
        BlockIDButton <- gbutton("Find Blocks", cont=BlockDetectWindow, expand = TRUE, handler=function(h,...){
          HybRIDS_Sessions[[sel]]$BlockDetectionParams$ManualThresholds <<- as.numeric(unlist(strsplit(svalue(ManualThresholdsEntry), ",")))
          HybRIDS_Sessions[[sel]]$BlockDetectionParams$AutoThresholds <<- svalue(AutoDetectCheck)
          HybRIDS_Sessions[[sel]]$BlockDetectionParams$ManualFallback <<- svalue(ManualFallbackCheck)
          HybRIDS_Sessions[[sel]]$BlockDetectionParams$SDstringency <<- as.numeric(svalue(SDstringencyGedit))
          if(svalue(FilterCheck)==TRUE){
            HybRIDS_Sessions[[sel]]$findBlocks(svalue(TripletFilter))
          } else {
            HybRIDS_Sessions[[sel]]$findBlocks()
          }
          dispose(BlockDetectWindow)
        })
        visible(BlockDetectWindow) <- TRUE
      },
    TestDateHandler =
      function(){
        sel <- svalue(AnalysisSwitcher)
        DatingWindow <- gwindow(paste("Block Significance testing and Dating for HybRIDS Session:", sel), visible=FALSE)
        MutRateLabel <- glabel("Mutation Rate:", cont=DatingWindow)
        MutGedit <- gedit(as.character(HybRIDS_Sessions[[sel]]$BlockDatingParams$MutationRate), width=10, cont=DatingWindow)
        PThreshLabel <- glabel("P-Value Threshold:", cont=DatingWindow)
        PTgedit <- gedit(as.character(HybRIDS_Sessions[[sel]]$BlockDatingParams$PValue), cont=DatingWindow, width=10)
        AllDateCheck <- gcheckbox("Apply Triplet Selecton Filter?", cont=DatingWindow)
        BonfCheck <- gcheckbox("Use Bonferroni correction?", cont=DatingWindow)
        DateBlockButton <- gbutton("Date Blocks", cont = DatingWindow, handler=function(h,...){
          HybRIDS_Sessions[[sel]]$BlockDatingParams$MutationRate <<- as.numeric(svalue(MutGedit))
          HybRIDS_Sessions[[sel]]$BlockDatingParams$PValue <<- as.numeric(svalue(PTgedit))
          if(svalue(AllDateCheck)==TRUE){
            HybRIDS_Sessions[[sel]]$dateBlocks(svalue(TripletFilter))
          } else {
            HybRIDS_Sessions[[sel]]$dateBlocks()
          }
          dispose(DatingWindow)
        })
        visible(DatingWindow) <- TRUE
      },
    CheckPresent =
      function(s){
        if(!s %in% names(HybRIDS_Sessions)){
          gmessage("A HybRIDS session of that name does not exist, perhaps you deleted it earlier but still have it typed in the analysis selector box", icon="error")
          stop("A HybRIDS session of that name does not exist")
        }
      },
    PlotExportHandler =
      function(){
        sel <- svalue(AnalysisSwitcher)
        CheckPresent(sel)
        DataWindow <- gwindow(paste("Data & Plotting of HybRIDS session:", sel), visible=FALSE)
        MainGroup <- ggroup(horizontal=TRUE, cont=DataWindow)
        PlotFrame <- gframe("Plotting Sequence Similarity", horizontal=TRUE, cont=MainGroup)
        Sub1 <- ggroup(horizontal=F, cont=PlotFrame)
        Sub2 <- ggroup(horizontal=F, cont=PlotFrame)
        barsbool <- "Bars" %in% HybRIDS_Sessions[[sel]]$PlottingParams$What
        linesbool <- "Lines" %in% HybRIDS_Sessions[[sel]]$PlottingParams$What
        WhatCheck1 <- gcheckbox("Plot Bars", checked=barsbool, cont=Sub1)
        WhatCheck2 <- gcheckbox("Plot Lines", checked=linesbool, cont=Sub1)
        TitleCheck <- gcheckbox("Add Title to Plot?", checked=HybRIDS_Sessions[[sel]]$PlottingParams$PlotTitle, cont=Sub1)
        CombinedTitleCheck <- gcheckbox("Main title for combined plots?", checked=HybRIDS_Sessions[[sel]]$PlottingParams$CombinedTitle, cont=Sub1)
        TitleSizeLabel <- glabel("Title Size", cont=Sub1)
        TitleSizeGedit <- gedit(as.character(HybRIDS_Sessions[[sel]]$PlottingParams$TitleSize), cont=Sub1)
        TitleFaceLabel <- glabel("Title Face", cont=Sub1)
        TitleFaceGedit <- gedit(as.character(HybRIDS_Sessions[[sel]]$PlottingParams$TitleFace), cont=Sub1)
        TitleColourLabel <- glabel("Title Colour:", cont=Sub1)
        TitleColourGedit <- gedit(as.character(HybRIDS_Sessions[[sel]]$PlottingParams$TitleColour), cont=Sub1)
        XFrame <- gframe("X Axis Settings", horizontal=F, cont=Sub2)
        XTitleCheck <- gcheckbox("Add Title to X axis?", checked=HybRIDS_Sessions[[sel]]$PlottingParams$XTitle, cont=XFrame)
        XTitleSizeLabel <- glabel("X Axis title size:", cont=XFrame)
        XTitleSizeGedit <- gedit(as.character(HybRIDS_Sessions[[sel]]$PlottingParams$XTitleFontSize), cont=XFrame)
        XTitleColourLabel <- glabel("X Axis title colour:", cont=XFrame)
        XTitleColourGedit <- gedit(HybRIDS_Sessions[[sel]]$PlottingParams$XTitleColour, cont=XFrame)
        XlabelsCheck <- gcheckbox("Add Labels to X axis?", checked=HybRIDS_Sessions[[sel]]$PlottingParams$XLabels, cont=XFrame)
        XLabelSizeLabel <- glabel("X Axis labels size", cont=XFrame)
        XLabelSizeGedit <- gedit(as.character(HybRIDS_Sessions[[sel]]$PlottingParams$XLabelSize), cont=XFrame)
        XLabelColourLabel <- glabel("X axis labels colour:", cont=XFrame)
        XLabelColourGedit <- gedit(HybRIDS_Sessions[[sel]]$PlottingParams$XLabelColour, cont=XFrame)
        YFrame <- gframe("Y Axis Labels", horizontal=FALSE, cont=Sub2)
        YTitleCheck <- gcheckbox("Add Title to Y axis?", checked=HybRIDS_Sessions[[sel]]$PlottingParams$YTitle, cont=YFrame)
        YTitleSizeLabel <- glabel("Y Axis title size:", cont=YFrame)
        YTitleSizeGedit <- gedit(as.character(HybRIDS_Sessions[[sel]]$PlottingParams$YTitleFontSize), cont=YFrame)
        YTitleColourLabel <- glabel("Y Axis title colour:", cont=YFrame)
        YTitleColourGedit <- gedit(HybRIDS_Sessions[[sel]]$PlottingParams$YTitleColour, cont=YFrame)
        YlabelsCheck <- gcheckbox("Add Labels to Y axis?", checked=HybRIDS_Sessions[[sel]]$PlottingParams$YLabels, cont=YFrame)
        YLabelSizeLabel <- glabel("Y Axis labels size", cont=YFrame)
        YLabelSizeGedit <- gedit(as.character(HybRIDS_Sessions[[sel]]$PlottingParams$YLabelSize), cont=YFrame)
        YLabelColourLabel <- glabel("Y axis labels colour:", cont=YFrame)
        YLabelColourGedit <- gedit(as.character(HybRIDS_Sessions[[sel]]$PlottingParams$YLabelColour), cont=YFrame)
        LegendCheck <- gcheckbox("Include Legend?", checked=HybRIDS_Sessions[[sel]]$PlottingParams$Legends, cont=Sub1)
        LegendFontSizeLabel <- glabel("Legend font size:", cont=Sub1)
        LegendFontSizeGedit <- gedit(as.character(HybRIDS_Sessions[[sel]]$PlottingParams$LegendFontSize), cont=Sub1)
        MosaicScaleLabel <- glabel("Mosaic Scale:", cont=Sub1)
        MosaicScaleGedit <- gedit(as.character(HybRIDS_Sessions[[sel]]$PlottingParams$MosaicScale), cont=Sub1)
        PlotButton <- gbutton("Create & View Plot", cont=Sub1, expand=TRUE, handler=function(h,...){
          if(svalue(WhatCheck1) == TRUE && svalue(WhatCheck2) == FALSE){
            HybRIDS_Sessions[[sel]]$PlottingParams$What <<- "Bars"
          }
          if(svalue(WhatCheck1) == FALSE && svalue(WhatCheck2) == TRUE){
            HybRIDS_Sessions[[sel]]$PlottingParams$What <<- "Lines"
          }
          if(svalue(WhatCheck1) == TRUE && svalue(WhatCheck2) == TRUE){
            HybRIDS_Sessions[[sel]]$PlottingParams$What <<- c("Bars", "Lines")
          }
          HybRIDS_Sessions[[sel]]$setParameters("Plotting", PlotTitle = svalue(TitleCheck), CombinedTitle = svalue(CombinedTitleCheck),
                                                TitleSize = as.numeric(svalue(TitleSizeGedit)), TitleFace = svalue(TitleFaceGedit),
                                                TitleColour = svalue(TitleColourGedit), XLabels = svalue(XlabelsCheck),
                                                YLabels = svalue(YlabelsCheck), XTitle = svalue(XTitleCheck),
                                                XTitleFontSize = as.numeric(svalue(XTitleSizeGedit)), XTitleColour = svalue(XTitleColourGedit),
                                                XLabelSize = as.numeric(svalue(XLabelSizeGedit)),
                                                XLabelColour = svalue(XLabelColourGedit), YTitle = svalue(YTitleCheck),
                                                YTitleFontSize = as.numeric(svalue(YTitleSizeGedit)), YTitleColour = svalue(YTitleColourGedit),
                                                YLabelSize = as.numeric(svalue(YLabelSizeGedit)), YLabelColour = svalue(YLabelColourGedit),
                                                Legends = svalue(LegendCheck), LegendFontSize = as.numeric(svalue(LegendFontSizeGedit)),
                                                MosaicScale = as.numeric(svalue(MosaicScaleGedit)))
          if(svalue(TripletFilter) == ""){
            gmessage("You need to specify a triplet to plot in the Triplet Selection Filter on the HybRIDS control panel e.g. Seq1:Seq2:Seq2", icon="error")
          } else {
            CurrentPlot <<- HybRIDS_Sessions[[sel]]$plotSS(svalue(TripletFilter))
            print(CurrentPlot)
          }
        })
        SavePlotButton <- gbutton("Save Current Plot", expand=TRUE, cont=Sub1, handler=function(h,...){
          gmessage("Choose a Filename for example MyPlot.png or MyPlot.pdf, The format will be detected from the filename you give e.g. '.pdf' = PDF format, or '.png' = PNG format", icon="info")
          fileName <- gfile(text="Choose a filename", 
                            type="save",
                            action="print",
                            handler= 
                              function(h,...){
                                do.call(h$action, list(h$file))
                              }
          )
          PlotSaveWin <- gwindow("Select Options for Plot Save", visible=FALSE)
          HeightLabel <- glabel("Image Height (inches):", cont=PlotSaveWin)
          HeightGedit <- gedit("10", cont=PlotSaveWin)
          WidthLabel <- glabel("Image Width (inches):", cont=PlotSaveWin)
          WidthGedit <- gedit("20", cont=PlotSaveWin)
          SaveButton <- gbutton("Save", cont=PlotSaveWin, expand=T, handler=function(h,...){
            ggsave(filename = fileName, plot = CurrentPlot, height=as.numeric(svalue(HeightGedit)), width=as.numeric(svalue(WidthGedit)))
            dispose(PlotSaveWin)
          })
          visible(PlotSaveWin) <- TRUE
        })
        BlocksFrame <- gframe("Detected Blocks", cont=MainGroup)
        BlocksButton <- gbutton("Display Table", cont=BlocksFrame)
        visible(DataWindow) <- TRUE
      }
  )
)