# Operations to be done upon package setup - mostly the GUI...


.onLoad <- function(...){
  
  options(guiToolkit="tcltk")
  # Define environments for the GUI guided session.
  HybRIDSenv <- new.env()
  HybRIDSenv$fastaseqs <- list()
  HybRIDSenv$analysisSessions <- list()
  
  # Define internal GUI functions.
  ## Internal function that selects the sequence set from the collection according to the choice made in the GUI.
  seqselector <- function(sequenceCollection, selectionChoice){
    output <- sequenceCollection[[which(names(sequenceCollection) == selectionChoice)]]
    return(output)
  }
  seqinputsorter <- function(newset, newsetname){
    HybRIDSenv$fastaseqs[[length(HybRIDSenv$fastaseqs)+1]] <- newset
    names(HybRIDSenv$fastaseqs)[[length(HybRIDSenv$fastaseqs)]] <- newsetname
    FastaSeqSelect[] <- names(HybRIDSenv$fastaseqs)
    names(sequencestodo) <- "Sequences"
  }
  deleteseqset <- function(valuetodel){
    svalue(FastaSeqSelect) <- ""
    HybRIDSenv$fastaseqs <- HybRIDSenv$fastaseqs[which(names(HybRIDSenv$fastaseqs) != valuetodel)]
    FastaSeqSelect[] <- FastaSeqSelect[FastaSeqSelect[] != valuetodel]
  }
  analysisselector <- function(analysissessions, selection){
    return(analysissessions[[which(names(analysissessions) == selection)]])
  }
  analysisindexer <- function(analysissessions, selection){
    return(which(names(analysissessions) == selection))
  }
  ssanalysissorter <- function(analysis, analysisname, sequencedep){
    if("HybRIDSseqsim" %in% class(analysis)){
      outputanalysis <- list(SSAnalysis=analysis, Blocks=as.HybRIDSblock(list()), SeqDep = sequencedep)
    } else {
      outputanalysis <- list(SSAnalysis=analysis, Blocks=as.HybRIDSblockSET(list(length(analysis))), SeqDep = sequencedep)
    }
    outputanalysis <- list(SSAnalysis=analysis, Blocks=list(length(analysis)), SeqDep = sequencedep)
    HybRIDSenv$analysisSessions[[length(HybRIDSenv$analysisSessions)+1]] <- outputanalysis
    names(HybRIDSenv$analysisSessions)[[length(HybRIDSenv$analysisSessions)]] <- analysisname
    AnalysisSwitcher[] <- names(HybRIDSenv$analysisSessions)
  }
  
  dothewindow <- function(){
    if("HybRIDSseqsim" %in% class(HybRIDSenv$analysisSessions[[analysisindexer(HybRIDSenv$analysisSessions, svalue(AnalysisSwitcher))]]$SSAnalysis)){
      svalue(WindowSizeLab) <- paste("Sliding Window Size: ", HybRIDSenv$analysisSessions[[analysisindexer(HybRIDSenv$analysisSessions, svalue(AnalysisSwitcher))]]$SSAnalysis[[1]]$WindowEnd[1] - HybRIDSenv$analysisSessions[[analysisindexer(HybRIDSenv$analysisSessions, svalue(AnalysisSwitcher))]]$SSAnalysis[[1]]$WindowStart[1])
    } else {
      if(svalue(Triplets3) != "Pair" && svalue(Triplets2) != "Individual" && length(unique(c(svalue(Triplets1),svalue(Triplets2),svalue(Triplets3)))) == 3){
        svalue(WindowSizeLab) <- paste("Sliding Window Size: ", HybRIDSenv$analysisSessions[[analysisindexer(HybRIDSenv$analysisSessions, svalue(AnalysisSwitcher))]]$SSAnalysis[[TripletIndexes(HybRIDSenv$analysisSessions[[analysisindexer(HybRIDSenv$analysisSessions, svalue(AnalysisSwitcher))]]$SSAnalysis, svalue(Triplets1), svalue(Triplets2), svalue(Triplets3))]][[1]]$WindowEnd[1] - HybRIDSenv$analysisSessions[[analysisindexer(HybRIDSenv$analysisSessions, svalue(AnalysisSwitcher))]]$SSAnalysis[[TripletIndexes(HybRIDSenv$analysisSessions[[analysisindexer(HybRIDSenv$analysisSessions, svalue(AnalysisSwitcher))]]$SSAnalysis, svalue(Triplets1), svalue(Triplets2), svalue(Triplets3))]][[1]]$WindowStart[1])
      }
    }
    if("HybRIDSseqsim" %in% class(HybRIDSenv$analysisSessions[[analysisindexer(HybRIDSenv$analysisSessions, svalue(AnalysisSwitcher))]]$SSAnalysis)){
      svalue(StepSizeLab) <- paste("Sliding Window Step Size: ", HybRIDSenv$analysisSessions[[analysisindexer(HybRIDSenv$analysisSessions, svalue(AnalysisSwitcher))]]$SSAnalysis[[1]]$WindowStart[2] - HybRIDSenv$analysisSessions[[analysisindexer(HybRIDSenv$analysisSessions, svalue(AnalysisSwitcher))]]$SSAnalysis[[1]]$WindowStart[1])
    } else {
      if(svalue(Triplets3) != "Pair" && svalue(Triplets2) != "Individual" && length(unique(c(svalue(Triplets1),svalue(Triplets2),svalue(Triplets3)))) == 3){
        svalue(StepSizeLab) <- paste("Sliding Window Step Size: ", HybRIDSenv$analysisSessions[[analysisindexer(HybRIDSenv$analysisSessions, svalue(AnalysisSwitcher))]]$SSAnalysis[[TripletIndexes(HybRIDSenv$analysisSessions[[analysisindexer(HybRIDSenv$analysisSessions, svalue(AnalysisSwitcher))]]$SSAnalysis, svalue(Triplets1), svalue(Triplets2), svalue(Triplets3))]][[1]]$WindowStart[2] - HybRIDSenv$analysisSessions[[analysisindexer(HybRIDSenv$analysisSessions, svalue(AnalysisSwitcher))]]$SSAnalysis[[TripletIndexes(HybRIDSenv$analysisSessions[[analysisindexer(HybRIDSenv$analysisSessions, svalue(AnalysisSwitcher))]]$SSAnalysis, svalue(Triplets1), svalue(Triplets2), svalue(Triplets3))]][[1]]$WindowStart[1])
      }
    }
  }
  
  # GUI Construction Starts now...
  MainWindow <- gwindow("HybRIDS: Simple & quick detection & dating of recombinant regions in DNA sequences", visible=FALSE)
  MainWinGroup <- ggroup(horizontal=TRUE, container=MainWindow)
  Notebook <- gnotebook(tab.pos = 3, closebuttons = FALSE, container = MainWinGroup, toolkit = guiToolkit())
  SecondMainwinGroup <- ggroup(horizontal=FALSE, container=MainWinGroup)
  DataSelection <- gframe("Data Selector", cont=SecondMainwinGroup, horizontal=F)
  DataSelectGroup2 <- gframe("\nAnalysis Selection", cont=DataSelection, horizontal=T)
  AnalysisSwitcher <- gcombobox("", cont=DataSelectGroup2, width=100, handler=function(h,...){
    if("HybRIDSseqsim" %in% class(analysisselector(HybRIDSenv$analysisSessions, svalue(AnalysisSwitcher))$SSAnalysis)){
      names <- analysisselector(HybRIDSenv$analysisSessions, svalue(AnalysisSwitcher))$SSAnalysis$ContigNames
    } else {
      if("HybRIDSseqsimSET" %in% class(analysisselector(HybRIDSenv$analysisSessions, svalue(AnalysisSwitcher))$SSAnalysis)){
        names <- unique(unlist(lapply(analysisselector(HybRIDSenv$analysisSessions, svalue(AnalysisSwitcher))$SSAnalysis, function(x) x$ContigNames)))
      }
    }
    Triplets1[] <- names
    Triplets2[] <- c("Individual", names)
    Triplets3[] <- c("Pair", names)
    if("HybRIDSseqsim" %in% class(HybRIDSenv$analysisSessions[[analysisindexer(HybRIDSenv$analysisSessions, svalue(AnalysisSwitcher))]]$SSAnalysis)){
      svalue(TripletNumLab) <- "Number of Triplets: 1"
    } else {
      svalue(TripletNumLab) <- paste("Number of Triplets: ", length(HybRIDSenv$analysisSessions[[analysisindexer(HybRIDSenv$analysisSessions, svalue(AnalysisSwitcher))]]$SSAnalysis))
    }
    dothewindow()
  })
  
  
  
  SwitcherButton <- gbutton("Delete Analysis", cont=DataSelectGroup2, handler=function(h,...){
    valuetodel <- svalue(AnalysisSwitcher)
    svalue(AnalysisSwitcher) <- ""
    HybRIDSenv$analysisSessions <- HybRIDSenv$analysisSessions[which(names(HybRIDSenv$analysisSessions) != valuetodel)]
    AnalysisSwitcher[] <- FastaSeqSelect[FastaSeqSelect[] != valuetodel]
  })
  # Data selection set of objects...
  DataSelectGroup <- gframe("Triplet Selection", cont=DataSelection, horizontal=T)
  Tg1 <- ggroup(cont=DataSelectGroup, horizontal=F)
  Tg2 <- ggroup(cont=DataSelectGroup, horizontal=F)
  Tg3 <- ggroup(cont=DataSelectGroup, horizontal=F)
  TLab1 <- glabel("First Sequence", cont=Tg1)
  TLab2 <- glabel("Second Sequence", cont=Tg2)
  TLab3 <- glabel("Third Sequence", cont=Tg3)
  Triplets1 <- gcombobox("", cont=Tg1, width=80, handler=function(h,...){dothewindow()})
  Triplets2 <- gcombobox("", cont=Tg2, width=80, handler=function(h,...){dothewindow()})
  Triplets3 <- gcombobox("", cont=Tg3, width=80, handler=function(h,...){dothewindow()})
  SelectedLabel <- glabel("The current triplet selection is: ", cont=DataSelection)
  DataSummary <- gframe("\nData Summary", cont=SecondMainwinGroup, expand=T)
  SequencesFrame <- gframe("Sequences: ", cont=DataSummary, horizontal=F)
  addSpring(DataSummary)
  vertanalysisframe <- ggroup(horizontal=F, cont=DataSummary)
  AnalysisFrame <- gframe("Analysis Session:", cont=vertanalysisframe, horizontal=F, expand=T)
  TripletFrame <- gframe("Individual Triplet:", cont=vertanalysisframe, horizontal=F, expand=T)
  SeqNameLabel <- glabel(text="Set Name: No Set Selected", cont=SequencesFrame)
  AlLengthLabel <- glabel(text="Alignment Length: 0 (0)", cont=SequencesFrame)
  SequenceNumLabel <- glabel(text="Number of Sequences: 0", cont=SequencesFrame)
  TripletTotalLabel <- glabel(text="Max Number of Triplets: 0", cont=SequencesFrame)
  TripletNumLab <- glabel("Number of Triplets: 0", cont=AnalysisFrame)
  WindowSizeLab <- glabel("Sliding Window Size: 0", cont=AnalysisFrame)
  StepSizeLab <- glabel("Step Size: 0", cont=AnalysisFrame)
  
  BlockDisplay <- gbutton("Show Blocks", cont = AnalysisFrame, handler = function(h,...){
    if(length(HybRIDSenv$analysisSessions) < 1){
      gmessage("You don't have any analyses saved in the HybRIDS GUI environment",icon="error")
      stop("You don't have any analyses saved in the HybRIDS GUI environment")
    }
    analysisindex <- analysisindexer(HybRIDSenv$analysisSessions, svalue(AnalysisSwitcher))
    if("HybRIDSblock" %in% class(HybRIDSenv$analysisSessions[[analysisindex]]$Blocks) || "HybRIDSdatedBlocks" %in% class(HybRIDSenv$analysisSessions[[analysisindex]]$Blocks)){
      BlocksTable <- SimplifyBlocks(HybRIDSenv$analysisSessions[[analysisindex]]$Blocks)
      blockwin <- gwindow(title=paste("Blocks for ", svalue(Triplets1), svalue(Triplets2), svalue(Triplets3)), visible=F)
      gtable(items=BlocksTable, container=blockwin)
      buttonsgroup <- gframe("Options", cont=blockwin, horizontal=T)
      savecsv <- gbutton("Save as .CSV", cont=buttonsgroup, handler=function(h,...){
        filename <- gfile(text="Choose a file", 
                          type="save",
                          filter = list("Fasta Files" = list(patterns = c("*.fasta","*.fas"))),
                          action="print",
                          handler= 
                            function(h,...){
                              do.call(h$action, list(h$file))
                            }
        )
        write.csv(BlocksTable, file=filename)
      })
      visible(blockwin) <- T
    } else {
      if(length(unique(c(svalue(Triplets1), svalue(Triplets2), svalue(Triplets3)))) > 2){
        BlocksTable <- SimplifyBlocks(HybRIDSenv$analysisSessions[[analysisindex]]$Blocks, svalue(Triplets1), svalue(Triplets2), svalue(Triplets3), returnEach=F, selectionOnly=T)
        blockwin <- gwindow(title=paste("Blocks for ", svalue(Triplets1), svalue(Triplets2), svalue(Triplets3)), visible=F)
        gtable(items=BlocksTable, container=blockwin)
        buttonsgroup <- gframe("Options", cont=blockwin, horizontal=T)
        savecsv <- gbutton("Save as .CSV", cont=buttonsgroup, handler=function(h,...){
          filename <- gfile(text="Choose a file", 
                            type="save",
                            filter = list("Fasta Files" = list(patterns = c("*.fasta","*.fas"))),
                            action="print",
                            handler= 
                              function(h,...){
                                do.call(h$action, list(h$file))
                              }
          )
          write.csv(BlocksTable, file=filename)
        })
        visible(blockwin) <- T
      } else {
        gmessage("Enter a valid triplet selection", title = "Error: Invalid Triplet Selection")
      }
    }
  })
  addHandlerMouseMotion(BlockDisplay, handler = function(h,...){
    svalue(HintBar) <- "Display all blocks detected for the triplet/pair selected for in the Data Select pane as a dataframe/table. Appears in new window.."
  })
  
  AboutButton <- gbutton("About", cont=SecondMainwinGroup, handler=function(h,...){
    Aboutwin <- gwindow("About HybRIDS", visible=F)
    aboutgroup1 <- ggroup(cont = Aboutwin)
    aboutgroup2 <- ggroup(cont=aboutgroup1, horizontal=F)
    logo <- gimage(filename = "HybRIDSlogo.png", dirname = "inst", container = aboutgroup2, expand=T)
    glabel("HybRIDS - recombination detection, dating &\n plotting of mosaic genome structures.\n", cont=aboutgroup2)
    glabel("Ben J. Ward, Cock van Oosterhout, Mark McMullan\n", cont=aboutgroup2)
    glabel("University of East Anglia\nSchool of Environmental Sciences\nELSA - Earth & Life Systems Alliance\n", cont=aboutgroup2)
    glabel("Citations:\n", cont=aboutgroup2)
    visible(Aboutwin) <- TRUE
  })
  # Write some code for the Helpful hints bar
  HintBar <- gstatusbar("Ready...", cont=MainWindow)
  
  # DNA Data preparation box of GUI
  framevert <- ggroup(cont=Notebook, horizonal=T, expand = TRUE, label="1). Sequence Preparation & Analysis")
  frame0 <- ggroup(cont=framevert, horizontal=F)
  subframe0 <- gframe("Read and Load DNA FASTA file", cont = frame0, horizontal = FALSE)
  preanalysis1 <- gcheckbox("Pre-Sort by Raw Sequence Similarity Threshold", checked = FALSE, use.togglebutton=FALSE, handler= NULL, action = NULL, container = subframe0)
  addHandlerMouseMotion(preanalysis1, handler = function(h,...){
    svalue(HintBar) <- "Select this option and set the threshold below it. Sequences with raw distances less than the threshold will not be able to form triplets."
  })
  
  labpre <- glabel("Pre-Sort Raw Sequence Similarity Threshold", cont=subframe0)
  prethreshold <- gedit(text="0.01", width = 20, cont=subframe0)
  preanalysis2 <- gcheckbox("Pre-Sort by Distance Distribution Dependent Cutoff", checked = FALSE, use.togglebutton=FALSE, handler= NULL, action = NULL, container = subframe0)
  addHandlerMouseMotion(preanalysis2, handler = function(h,...){
    svalue(HintBar) <- "Select this option to prevent sequences with small distances that fall outside the norm from forming triplets. HybRIDS will figure out what's outside the norm..."
  })
  DataButton <- gbutton("Select and Load DNA Data (FASTA)", container=subframe0, expand=TRUE,
                        handler=function(h,...){
                          sequencesname <- ginput(message="Give a name to the sequence set you're loading in", title="Sequence Names")
                          filename <- gfile(text="Choose a file", 
                                            type="open",
                                            filter = list("Fasta Files" = list(patterns = c("*.fasta","*.fas"))),
                                            action="print",
                                            handler= 
                                              function(h,...){
                                                do.call(h$action, list(h$file))
                                              }
                          )
                          if(svalue(preanalysis1) == TRUE && svalue(preanalysis2) == FALSE){
                            currentseqs <- dna.data.prepare(filename, 2, as.numeric(svalue(prethreshold)))
                            seqinputsorter(currentseqs, sequencesname)
                          } else {
                            if(svalue(preanalysis2) == TRUE && svalue(preanalysis1) == FALSE){
                              currentseqs <- dna.data.prepare(filename, 3, as.numeric(svalue(prethreshold)))
                              seqinputsorter(currentseqs, sequencesname)
                            } else {
                              if(svalue(preanalysis2) == FALSE && svalue(preanalysis1) == FALSE){
                                currentseqs <- dna.data.prepare(filename, 1, as.numeric(svalue(prethreshold)))
                                seqinputsorter(currentseqs, sequencesname)
                              } else {
                                gmessage("You need to select a valid method of pre-analysis - you can either use one of the methods or none. Not both.
                                         Unselect one.", title="Invalid Pre-Analysis Selection", icon="error")
                              }
                              }
                            }
                          svalue(HintBar) <- paste("DNA Data Loaded as ", sequencename, "...", sep="")
                          })
  addHandlerMouseMotion(DataButton, handler = function(h,...){
    svalue(HintBar) <- "Click to load in the sequences from fasta file..."
  })
  SequenceSelectGroup <- gframe("\nSelect Sequence Set", cont=subframe0, horizontal=T, expand=T)
  FastaSeqSelect <- gcombobox("", width = 100, cont=SequenceSelectGroup, editable=TRUE, handler=function(h,...){
    if(!svalue(FastaSeqSelect) %in% names(HybRIDSenv$fastaseqs)){
      sequencestodo[] <- ""
    } else {
      sequencestodo[] <- c("Auto", seqselector(HybRIDSenv$fastaseqs, svalue(FastaSeqSelect))$ContigNames)
    }
    names(sequencestodo) <- "Sequences"
    if(is.na(svalue(FastaSeqSelect)) || svalue(FastaSeqSelect) == ""){
      svalue(SeqNameLabel) <- "Set Name: "
      svalue(AlLengthLabel) <- "Alignment Length: "
      svalue(SequenceNumLabel) <- "Number of Sequences: "
      svalue(TripletTotalLabel) <- "Max Number of Triplets: "
    } else {
      svalue(SeqNameLabel) <- paste("Set Name: ", svalue(FastaSeqSelect))
      svalue(AlLengthLabel) <- paste("Alignment Length: ", ncol(seqselector(HybRIDSenv$fastaseqs, svalue(FastaSeqSelect))$Sequence))
      svalue(SequenceNumLabel) <- paste("Number of Sequences: ", nrow(seqselector(HybRIDSenv$fastaseqs, svalue(FastaSeqSelect))$Sequence))
      svalue(TripletTotalLabel) <- paste("Max Number of Triplets: ", length(combn(1:nrow(seqselector(HybRIDSenv$fastaseqs, svalue(FastaSeqSelect))$Sequence), 3, simplify=F)))
    }
  })
  addHandlerMouseMotion(FastaSeqSelect, handler = function(h,...){
    svalue(HintBar) <- "Select the currently active DNA sequence set..."
  })
  addSpring(SequenceSelectGroup)
  Fastabutton <- gbutton("Delete Sequence Set", cont=SequenceSelectGroup, handler=function(h,...){
    valuetodel <- svalue(FastaSeqSelect)
    if(any(valuetodel == unlist(lapply(HybRIDSenv$analysisSessions, function(x) x$SeqDep)))){
      delconf <- gconfirm("You have analysis sessions which will depend on these sequences,\nfor dating of identified blocks.\nAre yousure you want to delete the sequences?")
      if(delconf == TRUE){
        deleteseqset(valuetodel)
      }
    } else {
      deleteseqset(valuetodel)
    }
  })
  addHandlerMouseMotion(preanalysis2, handler = function(h,...){
    svalue(HintBar) <- "Delete the currently active sequence set..."
  })
  
  
  
  
  # Sequence Similarity Analysis Options
  subframe1 <- gframe("\nAnalysis Options", container=frame0)
  winoptframe <- ggroup(cont=subframe1, horizontal=F)
  subwinopt <- gframe("Window Settings", cont = winoptframe, horizontal=F)
  lab3 <- glabel("Window Size",cont=subwinopt)
  winsizespin <- gedit(text="100", width = 20, cont=subwinopt)
  addHandlerMouseMotion(winsizespin, handler = function(h,...){
    svalue(HintBar) <- "Set the size of the sliding window..."
  })
  lab4 <- glabel("Step Size",cont=subwinopt)
  stepsizespin <- gedit(text="1", width = 20, cont=subwinopt)
  addHandlerMouseMotion(stepsizespin, handler = function(h,...){
    svalue(HintBar) <- "Set the step size of the window..."
  })
  AnalysisTypeRadio <- gradio(c("New Analysis Session", "Re-Run Full Analysis", "Re-Run Analysis for Triplet"),
                              cont=winoptframe)
  addHandlerMouseMotion(AnalysisTypeRadio, handler = function(h,...){
    svalue(HintBar) <- "If Re-Run Analysis for Triplet is selected, the analysis of the triplet selected in the Data Selector pane will be re-run..."
  })
  
  Analyze1 <- gbutton("Analyse Sequence Similarity", container=winoptframe,
                      handler=function(h,...){
                        if(svalue(AnalysisTypeRadio) == "New Analysis Session"){
                          analysisname <- ginput(message="Give a name to the analysis you're creating", title="Name Analysis")
                          if(analysisname == "" || is.na(analysisname)){
                            gmessage("You need to give a valid name to the analysis!", title="Name your analysis!", icon="error")
                          } else {
                            if(all(!is.na(svalue(sequencestodo))) && length(svalue(sequencestodo)) == 1 && svalue(sequencestodo) == "Auto"){
                              analysis <- Calc.Similarity(seqselector(HybRIDSenv$fastaseqs, svalue(FastaSeqSelect)), window.size=as.numeric(svalue(winsizespin)), step.size=as.numeric(svalue(stepsizespin)))
                              ssanalysissorter(analysis, analysisname, svalue(FastaSeqSelect))
                            } else {
                              if(length(svalue(sequencestodo)) > 2 && !"Auto" %in% svalue(sequencestodo) && all(!is.na(svalue(sequencestodo)))){
                                analysisinput <- seqselector(HybRIDSenv$fastaseqs, svalue(FastaSeqSelect))
                                analysisinput <- SubSeq(analysisinput, svalue(sequencestodo))
                                analysis <- Calc.Similarity(analysisinput, window.size=as.numeric(svalue(winsizespin)), step.size=as.numeric(svalue(stepsizespin)))
                                ssanalysissorter(analysis, analysisname, svalue(FastaSeqSelect))
                              } else {
                                gmessage("You need to select a valid set of sequences from the Sequence list.\n\nSelect at least three sequences from the list or select 'All'.\n\n(Ctrl + click for multiple selection)", title="Invalid Sequence Selection", icon="error")
                              }
                            }
                          }
                        } else {
                          if(svalue(AnalysisTypeRadio) == "Re-Run Full Analysis"){
                            if(analysisselector(HybRIDSenv$analysisSessions, svalue(AnalysisSwitcher))$SeqDep %in% names(HybRIDSenv$fastaseqs)){
                              oldanalysis <- analysisselector(HybRIDSenv$analysisSessions, svalue(AnalysisSwitcher))
                              sequences <- seqselector(HybRIDSenv$fastaseqs, oldanalysis$SeqDep)
                              subseq <- SubSeq(sequences, oldanalysis$SSAnalysis$ContigNames)
                              analysis <- Calc.Similarity(subseq, window.size=as.numeric(svalue(winsizespin)), step.size=as.numeric(svalue(stepsizespin)))
                              HybRIDSenv$analysisSessions[[analysisindexer(HybRIDSenv$analysisSessions, svalue(AnalysisSwitcher))]] <- list(SSAnalysis = analysis, Blocks = list(), SeqDep = oldanalysis$SeqDep) 
                            } else {
                              gmessage("The FASTA sequences this analysis depends on cannot be found.", icon="error")
                            } 
                          } else {
                            if(svalue(AnalysisTypeRadio) == "Re-Run Analysis for Triplet"){
                              oldanalysisindex <- analysisindexer(HybRIDSenv$analysisSessions, svalue(AnalysisSwitcher))
                              sequencedep <- HybRIDSenv$analysisSessions[[oldanalysisindex]]$SeqDep
                              if(sequencedep %in% names(HybRIDSenv$fastaseqs)){
                                sequenceindex <- analysisindexer(HybRIDSenv$fastaseqs, sequencedep)
                                tripletchoice <- c(svalue(Triplets1), svalue(Triplets2), svalue(Triplets3))
                                tripletindexes <- TripletIndexes(HybRIDSenv$analysisSessions[[oldanalysisindex]]$SSAnalysis, sequence1=svalue(Triplets1), sequence2=svalue(Triplets2), sequence3=svalue(Triplets3))
                                oldsequence <- HybRIDSenv$fastaseqs[[sequenceindex]]
                                oldsequence$Combinations <- lapply(HybRIDSenv$analysisSessions[oldanalysisindex]$SSAnalysis[tripletindexes], function(x) which(HybRIDSenv$fastaseqs[[sequenceindex]]$ContigNames %in% x$ContigNames))
                                analyses <- Calc.Similarity(oldsequence, window.size=as.numeric(svalue(winsizespin)), step.size=as.numeric(svalue(stepsizespin)))
                                for(i in 1:length(oldsequence$Combinations)){
                                  HybRIDSenv$analysisSessions[[oldanalysisindex]]$SSAnalysis[[tripletindexes[i]]] <- analyses[[i]]
                                }
                              } else {
                                gmessage("You've deleted the sequence set that this analysis depends on from the workspace.\n\nYou cannot redo the SSAnalysis or identify and date blocks.\n\nYou can plot the analysis though.")
                              }
                            }
                          }
                        }
                        gmessage("Finished Analysing All Sequence Triplets", title="Finished Analysis", icon="info")
                        svalue(HintBar) <- "Finished Sequence Similarity Analysis..."
                      })
  addHandlerMouseMotion(Analyze1, handler = function(h,...){
    svalue(HintBar) <- "Run the analysis with selected options..."
  })
  
  
  
  sequencestodo <- gtable("", multiple=TRUE, cont=subframe1, expand=FALSE, width=100)
  addHandlerMouseMotion(sequencestodo, handler = function(h,...){
    svalue(HintBar) <- "Select which sequences to use in the analysis, use ctrl(cmd)+click for multiple selections. 'Auto' means all triplets not excluded by any pre-analysis done when loading in the fasta file..."
  })
  names(sequencestodo) <- "Sequences"
  
  
  
  
  
  # Block Detection Group
  detdateframe <- gframe("Block Detection & Dating", cont=framevert, horizontal=F)
  blockdetframe <- gframe("\nBlock Detection", cont=detdateframe)
  threshframe <- ggroup(cont=blockdetframe, horizontal=F)
  thresh2frame <- ggroup(cont=blockdetframe, horizontal=F)
  mutratelab <- glabel("Manual Similarity Thresholds", cont=threshframe)
  threshframe2 <- ggroup(cont=threshframe, horizontal=T, expand=T)
  addman <- gedit(width=5, cont=threshframe2)
  addHandlerMouseMotion(addman, handler = function(h,...){
    svalue(HintBar) <- "Specify new manual thresholds by entering them in here and pressing the add button..."
  })
  perclab <- glabel("%", cont=threshframe2)
  butgroup <- ggroup(cont=threshframe2, horizontal=F)
  addbut <- gbutton("Add", cont=butgroup, handler=function(h,...){
    if(!svalue(addman) %in% manvals[]){
      curvals <- manvals[]
      newvals <- svalue(addman)
      manvals[] <- c(curvals,newvals)
      names(manvals) <- "Man Thresholds"
    } else {
      gmessage("You already have that manual threshold in the list!", icon="error")
    }
  })
  delbut <- gbutton("Delete", cont=butgroup, handler=function(h,...){
    manvals[] <- manvals[manvals[] != svalue(manvals)]
    names(manvals) <- "Man Thresholds"
  })
  manvals <- gtable(90, multiple=TRUE, width = 100, height=100, cont=thresh2frame)
  addHandlerMouseMotion(manvals, handler = function(h,...){
    svalue(HintBar) <- "To remove a manual value, click it and then the Delete button, to multiple select manual values for analysis, use ctrl(cmd)+click for multi-select..."
  })
  names(manvals) <- "Man Thresholds"
  autodetectcheck <- gcheckbox("Auto-detect Thresholds?", checked = TRUE, use.togglebutton=FALSE, handler = NULL, action = NULL,
                               container = threshframe)
  addHandlerMouseMotion(autodetectcheck, handler = function(h,...){
    svalue(HintBar) <- "Select to have HybRIDS autodetect thresholds and try to decide which regions are putative blocks..."
  })
  
  
  
  manfallcheck <- gcheckbox("Fallback to manual values?", checked = TRUE, use.togglebutton = FALSE, handler = NULL, action = NULL,
                            container = threshframe)
  addHandlerMouseMotion(manfallcheck, handler = function(h,...){
    svalue(HintBar) <- "If TRUE, then if HybRIDS fails to find any thresholds it finds suitable it can fallack to the manual values you specify - this is recommended generally..."
  })
  allcheck <- gcheckbox("Run for All Triplets?", checked = TRUE, cont=threshframe)
  addHandlerMouseMotion(allcheck, handler = function(h,...){
    svalue(HintBar) <- "Set to true to run for all triplets, otherwise only the one selected for by the Data Selector pane will be used..."
  })
  
  sdstringencylab <- glabel("SD Stringency", cont=threshframe)
  sdstringency <- gedit(text = 2, width=5, cont=threshframe)
  addHandlerMouseMotion(sdstringency, handler = function(h,...){
    svalue(HintBar) <- "With lower values, HybRIDS will be more strict deciding what thresholds may be valid. The default is to consider thresholds half the SD above the mean SS as potentially valid (value of 2)..."
  })
  
  blockidbutton <- gbutton("Find Blocks", cont=thresh2frame, handler=function(h,...){
    if(is.na(svalue(manvals))){
      gmessage("You have not selected a manual threshold to either use to detect blocks\nor to use as a fallback for automatic threshold detection, having at least one is advised,\neven if only for the fallback.", title="No Manual Threshold Set", icon="error")
    } else {
      analysisindex <- analysisindexer(HybRIDSenv$analysisSessions, svalue(AnalysisSwitcher))
      if(svalue(allcheck) == TRUE || "HybRIDSseqsim" %in% class(HybRIDSenv$analysisSessions[[analysisindex]]$SSanalysis)){
        Blocks <- Identify.Blocks(HybRIDSenv$analysisSessions[[analysisindex]]$SSAnalysis, manual.thresholds=as.numeric(svalue(manvals)), autodetectThresholds=svalue(autodetectcheck), use.manual.fallback=svalue(manfallcheck), standard.deviation.stringency=as.numeric(svalue(sdstringency)))
        HybRIDSenv$analysisSessions[[analysisindex]]$Blocks <- Blocks
      } else {
        tripletindicies <- TripletIndexes(HybRIDSenv$analysisSessions[[analysisindex]]$SSAnalysis, svalue(Triplets1), svalue(Triplets2), svalue(Triplets3))
        subsetanalysis <- HybRIDSenv$analysisSessions[[analysisindex]]$SSAnalysis[[tripletindicies]]
        Blocks <- Identify.Blocks(subsetanalysis, manual.thresholds=as.numeric(svalue(manvals)), autodetectThresholds=T, use.manual.fallback=T, standard.deviation.stringency=2)
        for(i in 1:length(analysis)){
          HybRIDSenv$analysisSessions[[analysisindex]]$Blocks[[i]] <- Blocks[[i]]
        }
      }
    }
    gmessage("Finished Block Detection")
    svalue(HintBar) <- "Finished Block Detection..."
  })
  addHandlerMouseMotion(blockidbutton, handler = function(h,...){
    svalue(HintBar) <- "Start finding potential blocks..."
  })
  
  # Section for Block Dating.
  DatingFrame <- gframe("\nDating Blocks", cont=detdateframe, expand=T, horizontal=FALSE)
  mrategroup <- ggroup(cont=DatingFrame, horizontal=T, expand=T)
  MutRateLabel <- glabel("Mutation Rate:", cont=mrategroup)
  Mutgedit <- gedit("10e-8", width=10, cont=mrategroup)
  addSpring(mrategroup)
  AllDateCheck <- gcheckbox("Do for all triplets", cont=DatingFrame)
  DateBlockButton <- gbutton("Date Blocks", cont = mrategroup, handler = function(h,...){
    if(length(HybRIDSenv$analysisSessions) < 1){
      gmessage("No analyses sessions found!", icon="error")
      stop("No analysis sessions found")
    } 
    analysisindex <- analysisindexer(HybRIDSenv$analysisSessions, svalue(AnalysisSwitcher))
    if(!HybRIDSenv$analysisSessions[[analysisindex]]$SeqDep %in% names(HybRIDSenv$fastaseqs)){
      stop("DNA Sequence needed does not exist... did you delete it?")
    } else {
      dnatouse <- HybRIDSenv$fastaseqs[[which(names(HybRIDSenv$fastaseqs)==HybRIDSenv$analysisSessions[[analysisindex]]$SeqDep)]]
      if(svalue(AllDateCheck) == TRUE){
        cat(class(dnatouse))
        cat(class(HybRIDSenv$analysisSessions[[analysisindex]]$Blocks))
        HybRIDSenv$analysisSessions[[analysisindex]]$Blocks <- Estimate.Ages(HybRIDSenv$analysisSessions[[analysisindex]]$Blocks, dnatouse, as.numeric(svalue(Mutgedit)))
        
      }
      
      
    }
  })
  
  addHandlerMouseMotion(Mutgedit, handler = function(h,...){
    svalue(HintBar) <- "Adjust mutation rate used in dating function, a higher mutation rate means the function considers accumulation of mutations occurs faster and this affects the date..."
  })
  addHandlerMouseMotion(DateBlockButton, handler = function(h,...){
    svalue(HintBar) <- "Date the blocks..."
  })
  addHandlerMouseMotion(AllDateCheck, handler = function(h,...){
    svalue(HintBar) <- "If true, all blocks for all triplets will be dated, otherwise the blocks for the triplet selected will be dated..."
  })
  svalue(AllDateCheck) <- TRUE
  
  
  # Plotting Page of the notebook in the GUI.
  plotgroup1 <- ggroup(cont=Notebook, horizontal=T, label="2). Plotting & Visualisation")
  plotgroup2 <- ggroup(cont=plotgroup1, horizontal=F)
  TripletPlotFrame <- gframe("Plot Individual Triplets:", cont=plotgroup2, horizontal=F)
  PlotSSAnalysisFrame <- gframe("Plot Sequence Similariy Analysis:", cont=TripletPlotFrame)
  CheckGroup <- ggroup(cont=PlotSSAnalysisFrame, horizontal=F)
  LegendCheck <- gcheckbox("Include Legend?", checked=TRUE, cont=CheckGroup)
  ThresholdDensityCheck <- gcheckbox("Density Plot?", checked=FALSE, cont=CheckGroup)
  addHandlerMouseMotion(ThresholdDensityCheck, handler = function(h,...){
    svalue(HintBar) <- "Select to plot the threshold density for a triplet..."
  })
  OneCanvasCheck <- gcheckbox("On one Canvas?", checked=TRUE, cont=CheckGroup)
  addHandlerMouseMotion(OneCanvasCheck, handler = function(h,...){
    svalue(HintBar) <- "If FALSE, multiple plot images will be produced..."
  })
  TPvertgroup <- gframe("Axis Labelling:",cont=PlotSSAnalysisFrame, horizontal=F)
  AxisbpRadio <- gradio(c("Annotate x axis by bp", "Annotate x axis by window number"), cont=TPvertgroup)
  BpFreqGroup <- ggroup(cont=TPvertgroup, horizontal=T)
  BpFreqLabel <- glabel("x Axis Base Pair Frequency:", cont=BpFreqGroup)
  addSpring(BpFreqGroup)
  BpFreqScale <- gedit(cont=BpFreqGroup, width=5)
  addHandlerMouseMotion(BpFreqScale, handler = function(h,...){
    svalue(HintBar) <- "Modify this value to affect the spacing of the x axis ticks to make the plot neater if too many ticks are ruining the image..."
  })
  addHandlerMouseMotion(BpFreqLabel, handler = function(h,...){
    svalue(HintBar) <- "Modify this value to affect the spacing of the x axis ticks to make the plot neater if too many ticks are ruining the image..."
  })
  MosaicVertGroup <- gframe("Mosaic Bar Options:", cont=PlotSSAnalysisFrame, horizontal=F, expand=T)
  addHandlerMouseMotion(MosaicVertGroup, handler = function(h,...){
    svalue(HintBar) <- "Mosaic Bars are produced by colour mixing based on SS values on chunks of sequence..."
  })
  MosaicBarsCheck <- gcheckbox("Plot Mosaic Bars", checked=TRUE, cont=MosaicVertGroup)
  CondenseBarsCheck <- gcheckbox("Condense Mosaic Bars\n(Recommended)", checked=TRUE, cont=MosaicVertGroup)
  addHandlerMouseMotion(CondenseBarsCheck, handler = function(h,...){
    svalue(HintBar) <- "If you don't have this checked then the allwindows will be plotted as a grid, if checked then groups of windows across a sequence are combined and plotted to form one bar..."
  })
  MosaicGroup <- ggroup(cont=MosaicVertGroup, horizontal=T, expand=T)
  MosaicLabel <- glabel("Mosaic Scale:", cont=MosaicGroup)
  addSpring(MosaicGroup)
  MosaicScale <- gedit(cont=MosaicGroup, width=5)
  svalue(MosaicScale) <- 50
  svalue(BpFreqScale) <- 500
  SSAnalysisPlotButton <- gbutton("Plot SS Analysis", cont=PlotSSAnalysisFrame, handler=function(h,...){
    analysisindex <- analysisindexer(HybRIDSenv$analysisSessions, svalue(AnalysisSwitcher))
    if("HybRIDSseqsim" %in% class(HybRIDSenv$analysisSessions[[analysisindex]]$SSAnalysis)){
      if(svalue(AxisbpRadio) == "Annotate x axis by bp"){
        plot(HybRIDSenv$analysisSessions[[analysisindex]]$SSAnalysis, linesplot=T, legends=svalue(LegendCheck), baseannotate=T, bpfreq=as.numeric(svalue(BpFreqScale)), densityplot=svalue(ThresholdDensityCheck), mosaic.bars=svalue(MosaicBarsCheck), mosaic.scale=as.numeric(svalue(MosaicScale)), condense.mosaics=svalue(CondenseBarsCheck), labfontsize=as.numeric(svalue(AxisFontSize)), legfontsize=as.numeric(svalue(LegendFontSize)), onecanvas=svalue(OneCanvasCheck))
      } else {
        if(svalue(AxisbpRadio) == "Annotate x axis by window number"){
          plot(HybRIDSenv$analysisSessions[[analysisindex]]$SSAnalysis, linesplot=T, legends=svalue(LegendCheck), baseannotate=F, bpfreq=as.numeric(svalue(BpFreqScale)), densityplot=svalue(ThresholdDensityCheck), mosaic.bars=svalue(MosaicBarsCheck), mosaic.scale=as.numeric(svalue(MosaicScale)), condense.mosaics=svalue(CondenseBarsCheck), labfontsize=as.numeric(svalue(AxisFontSize)), legfontsize=as.numeric(svalue(LegendFontSize)), onecanvas=svalue(OneCanvasCheck))
        }
      }
    } else {
      if(svalue(Triplets2) == "Individual"){
        gmessage("Invalid Triplet Selection!\n\nPlease select three different sequence names, or two followed by 'Pair'", title="Invalid Triplet Selection", icon="error")
      } else {
        if(svalue(Triplets2) != "Individual" && svalue(Triplets3) == "Pair"){
          plot(HybRIDSenv$analysisSessions[[analysisindex]]$SSAnalysis, svalue(Triplets1), svalue(Triplets2))
        } else {
          if("HybRIDSseqsimSET" %in% class(HybRIDSenv$analysisSessions[[analysisindex]]$SSAnalysis) && length(unique(c(svalue(Triplets1), svalue(Triplets2), svalue(Triplets3)))) == 3 && svalue(Triplets2) != "Individual" && svalue(Triplets3) != "Pair"){
            tripletindex <- TripletIndexes(HybRIDSenv$analysisSessions[[analysisindex]]$SSAnalysis, sequence1=svalue(Triplets1), sequence2=svalue(Triplets2), sequence3=svalue(Triplets3))
            if(svalue(AxisbpRadio) == "Annotate x axis by bp"){
              plot(HybRIDSenv$analysisSessions[[analysisindex]]$SSAnalysis[[tripletindex]], linesplot=T, legends=svalue(LegendCheck), baseannotate=T, bpfreq=as.numeric(svalue(BpFreqScale)), densityplot=svalue(ThresholdDensityCheck), mosaic.bars=svalue(MosaicBarsCheck), mosaic.scale=as.numeric(svalue(MosaicScale)), condense.mosaics=svalue(CondenseBarsCheck), labfontsize=as.numeric(svalue(AxisFontSize)), legfontsize=as.numeric(svalue(LegendFontSize)), onecanvas=svalue(OneCanvasCheck))
            } else {
              if(svalue(AxisbpRadio) == "Annotate x axis by window number"){
                plot(HybRIDSenv$analysisSessions[[analysisindex]]$SSAnalysis[[tripletindex]], linesplot=T, legends=svalue(LegendCheck), baseannotate=F, bpfreq=as.numeric(svalue(BpFreqScale)), densityplot=svalue(ThresholdDensityCheck), mosaic.bars=svalue(MosaicBarsCheck), mosaic.scale=as.numeric(svalue(MosaicScale)), condense.mosaics=svalue(CondenseBarsCheck), labfontsize=as.numeric(svalue(AxisFontSize)), legfontsize=as.numeric(svalue(LegendFontSize)), onecanvas=svalue(OneCanvasCheck))
              }
            }
          } 
        }
      }
    }
  })
  
  PlotBlocksAndDatesFrame <- gframe("Plot Blocks and Dates:", cont=TripletPlotFrame)
  glabel("TEST", cont=PlotBlocksAndDatesFrame)
  
  FontFrame <- gframe("Font Options:", cont=plotgroup2, horizontal=F)
  addHandlerMouseMotion(FontFrame, handler = function(h,...){
    svalue(HintBar) <- "Font options for axis and legends for all plots..."
  })
  FontGroup1 <- ggroup(cont=FontFrame, horizontal=T, expand=TRUE)
  LegendFontSizeLab <- glabel("Legend Font Size:", cont=FontGroup1)
  addSpring(FontGroup1)
  LegendFontSize <- gedit(width=5, cont=FontGroup1)
  FontGroup2 <- ggroup(cont=FontFrame, horizontal=T)
  AxisFontSizeLab <- glabel("Axis Label Font Size:", cont=FontGroup2)
  AxisFontSize <- gedit(width=5, cont=FontGroup2)
  svalue(LegendFontSize) <- svalue(AxisFontSize) <- 14
  
  
  
  
  
  
  
  visible(MainWindow) <- TRUE
  

}