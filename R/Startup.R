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
    #widgets
    MainWindow = "ANY",
    MainWinGroup = "ANY",
    Notebook = "ANY",
    # All the widgets for the sequence analysis frame.
    SeqAnalysisFrame = "ANY",
    frame0 = "ANY",
    subframe0 = "ANY",
    tripletsetframe = "ANY",
    tripletsetbutton = "ANY",
    preanalysis1 = "ANY",
    labpre = "ANY",
    prethreshold = "ANY",
    preanalysis2 = "ANY",
    DataButton = "ANY",
    subframe1 = "ANY",
    winoptframe = "ANY",
    subwinopt = "ANY",
    lab3 = "ANY",
    winsizespin = "ANY",
    lab4 = "ANY",
    stepsizespin = "ANY",
    TripletFilterCheck = "ANY",
    # Block Detection Group
    detdateframe = "ANY",
    blockdetframe = "ANY",
    threshframe = "ANY",
    threshframe2 = "ANY",
    addman = "ANY",
    mutratelab = "ANY",
    thresh2frame = "ANY",
    perclab = "ANY",
    butgroup = "ANY",
    addbut = "ANY",
    delbut = "ANY",
    manvals = "ANY",
    autodetectcheck = "ANY", 
    manfallcheck = "ANY",
    allcheck = "ANY",
    sdstringencylab = "ANY",
    sdstringency = "ANY",
    blockidbutton = "ANY",
    # Block Dating Stuff
    DatingFrame = "ANY",
    mrategroup = "ANY",
    ptgroup = "ANY",
    MutRateLabel = "ANY",
    Mutgedit = "ANY",
    PThreshLabel = "ANY",
    PTgedit = "ANY",
    AllDateCheck  = "ANY",
    DateBlockButton = "ANY",
    BonfCheck = "ANY",
    
    AnalysisSwitcher = "ANY",
    SecondMainwinGroup = "ANY",
    DataSelection = "ANY",
    DataSelectGroup2 = "ANY",
    SwitcherButton = "ANY",
    DataSelectGroup = "ANY",
    TripletFilter = "ANY",
    
    
    DataSummary = "ANY",
    
    
    
    
    SelectedLabel = "ANY",
    
    StatusBar = "ANY",
    # Data
    HybRIDS_Sessions = "list"
    
  ),
  methods = list(
    initialize = function()
    {
      ## Initializing all of the GUI components...
      MainWindow <<- gwindow("HybRIDS", visible = FALSE)
      MainWinGroup <<- ggroup(horizontal=TRUE, container=MainWindow)
      Notebook <<- gnotebook(tab.pos = 3, closebuttons = FALSE, container = MainWinGroup, toolkit = guiToolkit())
      SeqAnalysisFrame <<- ggroup(cont=Notebook, horizonal=T, expand = TRUE, label="1). Sequence Preparation & Analysis")
      frame0 <<- ggroup(cont=SeqAnalysisFrame, horizontal=F)
      subframe0 <<- gframe("Create a new HybRIDS session from sequence file", cont = frame0, horizontal = FALSE)
      DataButton <<- gbutton("Step 1). Select and Load DNA Data (FASTA)", container=subframe0, expand=TRUE, handler=function(h,...){LoadIn()})
      tripletsetframe <<- gframe("Triplet Settings", container=frame0, horizontal = FALSE)
      preanalysis1 <<- gcheckbox("Pre-Sort by Raw Sequence Similarity Threshold", checked = FALSE, use.togglebutton=FALSE, handler = function(h,...){
        if(svalue(preanalysis1)==TRUE){
          if(svalue(preanalysis2)==TRUE){
            svalue(preanalysis2) <<- FALSE
          }
          HybRIDS_Sessions[[svalue(AnalysisSwitcher)]]$TripletParams$Method <<- 2
        } else {
          HybRIDS_Sessions[[svalue(AnalysisSwitcher)]]$TripletParams$Method <<- 1
        }
      }, action = NULL, container = tripletsetframe)
      labpre <<- glabel("Pre-Sort Raw Sequence Similarity Threshold", cont=tripletsetframe)
      prethreshold <<- gedit(text="0.01", width = 20, cont=tripletsetframe, handler = function(h,...){
        HybRIDS_Sessions[[svalue(AnalysisSwitcher)]]$TripletParams$SortThreshold <<- as.numeric(svalue(prethreshold))
        })
      preanalysis2 <<- gcheckbox("Pre-Sort by Distance Distribution Dependent Cutoff", checked = FALSE, use.togglebutton=FALSE, handler = function(h,...){
        if(svalue(preanalysis2)==TRUE){
          if(svalue(preanalysis1)==TRUE){
            svalue(preanalysis1) <<- FALSE
          }
          HybRIDS_Sessions[[svalue(AnalysisSwitcher)]]$TripletParams$Method <<- 3
        } else {
          HybRIDS_Sessions[[svalue(AnalysisSwitcher)]]$TripletParams$Method <<- 1
        }
      }, action = NULL, container = tripletsetframe)
      tripletsetbutton <<- gbutton("Step 2). Generate Triplets", cont=tripletsetframe, handler=function(h,...){
        HybRIDS_Sessions[[svalue(AnalysisSwitcher)]]$makeTripletCombos()
        gmessage("Made all the triplets", icon="info")
      })
      subframe1 <<- gframe("Analysis Options", container=frame0)
      winoptframe <<- ggroup(cont=subframe1, horizontal=F)
      subwinopt <<- gframe("Window Settings", cont = winoptframe, horizontal=F)
      lab3 <<- glabel("Window Size",cont=subwinopt)
      winsizespin <<- gedit(text="100", width = 20, cont=subwinopt, handler=function(h,...){
        HybRIDS_Sessions[[svalue(AnalysisSwitcher)]]$SSAnalysisParams$WindowSize <<- as.numeric(svalue(winsizespin))
      })
      lab4 <<- glabel("Step Size",cont=subwinopt)
      stepsizespin <<- gedit(text="1", width = 20, cont=subwinopt, handler=function(h,...){
        HybRIDS_Sessions[[svalue(AnalysisSwitcher)]]$SSAnalysisParams$StepSize <<- as.numeric(svalue(stepsizespin))
      })
      TripletFilterCheck <<- gcheckbox("Apply Triplet Selection Filter",
                                  cont=winoptframe)
      Analyze1 <- gbutton("Step 3). Analyse Sequence Similarity", container=winoptframe, handler=function(h,...){
        if(svalue(TripletFilterCheck)==TRUE){
          HybRIDS_Sessions[[svalue(AnalysisSwitcher)]]$analyzeSS(svalue(TripletFilter))
        } else {
          HybRIDS_Sessions[[svalue(AnalysisSwitcher)]]$analyzeSS()
        }
        gmessage("Finished analyzing sequence similarities of Sequence Triplets", icon="info")
      })
      # Block Detection Stuff
      detdateframe <<- gframe("Block Detection & Dating", cont=SeqAnalysisFrame, horizontal=F)
      blockdetframe <<- gframe("\nBlock Detection", cont=detdateframe)
      threshframe <<- ggroup(cont=blockdetframe, horizontal=F)
      thresh2frame <<- ggroup(cont=blockdetframe, horizontal=F)
      mutratelab <<- glabel("Manual Similarity Thresholds", cont=threshframe)
      threshframe2 <<- ggroup(cont=threshframe, horizontal=T, expand=T)
      addman <<- gedit(width=5, cont=threshframe2)
      perclab <<- glabel("%", cont=threshframe2)
      butgroup <<- ggroup(cont=threshframe2, horizontal=F)
      addbut <<- gbutton("Add", cont=butgroup, handler=function(h,...){
        if(!svalue(addman) %in% manvals[]){
          curvals <- manvals[]
          newvals <- svalue(addman)
          manvals[] <<- c(curvals,newvals)
          names(manvals) <<- "Man Thresholds"
        } else {
          gmessage("You already have that manual threshold in the list!", icon="error")
        }
      })
      delbut <<- gbutton("Delete", cont=butgroup, handler=function(h,...){
        manvals[] <<- manvals[manvals[] != svalue(manvals)]
        names(manvals) <<- "Man Thresholds"
      })
      manvals <<- gtable(90, multiple=TRUE, width = 100, height=100, cont=thresh2frame)
      autodetectcheck <<- gcheckbox("Auto-detect Thresholds?", checked = TRUE, use.togglebutton = FALSE, handler = function(h,...){
        HybRIDS_Sessions[[svalue(AnalysisSwitcher)]]$BlockDetectionParams$AutoThresholds <<- svalue(autodetectcheck)
      }, action = NULL,
                                   container = threshframe)
      manfallcheck <<- gcheckbox("Fallback to manual values?", checked = TRUE, use.togglebutton = FALSE, handler = function(h,...){
        HybRIDS_Sessions[[svalue(AnalysisSwitcher)]]$BlockDetectionParams$ManualFallback <<- svalue(manfallcheck)
      }, action = NULL,
                                container = threshframe)
      allcheck <<- gcheckbox("Apply Triplet Filter", cont=threshframe)
      sdstringencylab <<- glabel("SD Stringency", cont=threshframe)
      sdstringency <<- gedit(text = 2, width=5, cont=threshframe, handler=function(h,...){
        HybRIDS_Sessions[[svalue(AnalysisSwitcher)]]$BlockDetectionParams$SDstringency <<- as.numeric(svalue(sdstringency))
      })
      blockidbutton <<- gbutton("Step 4). Find Blocks", cont=thresh2frame, handler=function(h,...){
        HybRIDS_Sessions[[svalue(AnalysisSwitcher)]]$BlockDetectionParams$ManualThresholds <<- as.numeric(svalue(manvals))
        if(svalue(allcheck)==TRUE){
          HybRIDS_Sessions[[svalue(AnalysisSwitcher)]]$findBlocks(svalue(TripletFilter))
        } else {
          HybRIDS_Sessions[[svalue(AnalysisSwitcher)]]$findBlocks()
        }
      })
      # Block Dating Stuff
      DatingFrame <<- gframe("\nDating Blocks & P-Value Calculation", cont=detdateframe, expand=T, horizontal=FALSE)
      mrategroup <<- ggroup(cont=DatingFrame, horizontal=T, expand=T)
      ptgroup <<- ggroup(cont=DatingFrame, horizontal=T, expand=T)
      MutRateLabel <<- glabel("Mutation Rate:", cont=mrategroup)
      Mutgedit <<- gedit("10e-8", width=10, cont=mrategroup)
      PThreshLabel <<- glabel("P-Value Threshold:", cont=ptgroup)
      PTgedit <<- gedit("0.05", cont=ptgroup, width=10)
      addSpring(mrategroup)
      AllDateCheck <<- gcheckbox("Apply Triplet Selecton Filter?", cont=DatingFrame)
      BonfCheck <<- gcheckbox("Use Bonferroni correction?", cont=DatingFrame)
      DateBlockButton <<- gbutton("Step 5). Date Blocks", cont = mrategroup)
     
      # Data Selection Stuff
      SecondMainwinGroup <<- ggroup(horizontal=FALSE, container=MainWinGroup)
      DataSelection <<- gframe("Data Selector", cont=SecondMainwinGroup, horizontal=F)
      DataSelectGroup2 <<- gframe("HybRIDS Session Selection", cont=DataSelection, horizontal=T)
      AnalysisSwitcher <<- gcombobox("", cont=DataSelectGroup2, width=100, editable=T)
      SwitcherButton <<- gbutton("Delete Analysis", cont=DataSelectGroup2, handler=function(h,...){DelSession(svalue(AnalysisSwitcher))})
      DataSelectGroup <<- gframe("Triplet Selection Filter", cont=DataSelection, horizontal=T)
      TripletFilter <<- gedit("ALL", width=20, cont=DataSelectGroup)
      
      DataSummary <<- gbutton("HybRIDS Session Summary", cont=SecondMainwinGroup, expand=T,
                              handler = function(h, ...){
                                AnalysisDisplay(svalue(AnalysisSwitcher))
                              })
      
      
      StatusBar <<- gstatusbar("Ready...", cont=MainWindow)
      
      
      
      visible(MainWindow) <<- TRUE
    },
    AnalysisDisplay =
      function(sel) {
        CheckPresent()
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
        HybRIDS_Sessions <<- append(HybRIDS_Sessions, HybRIDS$new(filename))
        names(HybRIDS_Sessions)[[length(HybRIDS_Sessions)]] <<- sequencesname
        AnalysisSwitcher[] <<- names(HybRIDS_Sessions)
        gmessage("Loaded DNA Data and created new HybRIDS session", icon="info")
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
    CheckPresent =
      function(){
        if (!svalue(AnalysisSwitcher) %in% names(HybRIDS_Sessions)){
          gmessage("No HybRIDS session by that name exists, most likely you deleted it but still have it typed in the HybRIDS session selection box.", title="No Session", icon="error")
          stop()
        }  
      }
    
  )
)



