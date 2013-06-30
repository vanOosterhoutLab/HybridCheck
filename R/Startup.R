# Operations to be done upon package setup.
# Last altered by Ben J. Ward on 02/04/2013.

.onLoad <- function(...){
  
  options(guiToolkit="tcltk")
  
  conf <- gconfirm(message = "Would you like the HybRIDS graphical interface?
           This will be useful to you if you are not used to coding in R.", title="Load the HybRIDS GUI?", icon = "question")
  
  if(conf==TRUE){
    
    HybRIDSenv <- new.env()
    
    dnaseqs <- "NOPE"
    
    MainWindow <- gwindow("HybRIDS: Simple & quick detection & dating of recombinant regions in DNA sequences", visible=FALSE)
    
    title <- glabel("This package is maintained by Ben J. Ward. <b.ward@uea.ac.uk>", cont = MainWindow)
    
    MainWinGroup <- ggroup(horizontal=TRUE, container=MainWindow)
    
    logo <- gimage(filename = "HybRIDSlogo.gif", dirname = "inst", container = MainWinGroup)
    
    SecondMainwinGroup <- ggroup(horizontal=FALSE, container=MainWinGroup)
    
    DataButton <- gbutton("Load DNA Data (FASTA)", container=MainWindow,
                          handler=function(h,...){
                            filename <- gfile(text="Choose a file", 
                                  type="open", 
                                  action="print",
                                  handler= 
                                    function(h,...){
                                      do.call(h$action, list(h$file))
                                    }
                            )
                            HybRIDSenv$dnaseqs <- dna.data.prepare(filename, svalue(preanalysis), as.numeric(svalue(prethreshold)))
                          })
    
    
    Analyze <- gbutton("Analyze!", container=MainWindow,
                       handler=function(h,...){
                         analysis <- analyze(HybRIDSenv$dnaseqs, mutation.rate=as.numeric(svalue(mutspin)), win.size=as.numeric(svalue(winsizespin)), step=as.numeric(svalue(stepsizespin)))
                         HybRIDSenv$analysis <- analysis
                       })
    
    
    Plotbutton <- gbutton("Plot!", container=MainWindow,
                          handler=function(h,...){
                            if(length(HybRIDSenv$analysis)==1){
                              plot(HybRIDSenv$analysis, linesplot=svalue(linecheck), lineplot.legend=svalue(linelegcheck), baseannotate=svalue(bpcheck), bpfreq=as.numeric(svalue(bpspin)), mosaic.bars=svalue(mosaiccheck), mosaic.scale=as.numeric(svalue(mosaicspin)), combine.plots=svalue(combinecheck), condense.mosaics=svalue(condensecheck), labfontsize=as.numeric(svalue(fontspin1)), legfontsize=as.numeric(svalue(fontspin2)))
                            }
                          })
    
    
    frame1 <- gframe("Block Detection & Dating", cont=SecondMainwinGroup, horizontal=F)
    
    preanalysis <- gcheckbox("Pre-Analysis", checked = FALSE, use.togglebutton=FALSE, handler= NULL, action = NULL, container = frame1)
    
    labpre <- glabel("Pre-Analysis Distance Threshold", cont=frame1)
    
    prethreshold <- gedit(text="1", width = 20, cont=frame1) 
    
    lab1 <- glabel("Approximate Contig Length", cont=frame1)
    
    conlengthspin <- gedit(text="30000", width = 20, cont=frame1)
    
    lab2 <- glabel("Mutation Rate", cont=frame1)
    
    mutspin <- gedit(text="10e-8", width = 20, cont=frame1)
    
    frame2 <- gframe("Sliding Window", cont=SecondMainwinGroup, horizontal=F)
    
    lab3 <- glabel("Window Size",cont=frame2)
    
    winsizespin <- gedit(text="100", width = 20, cont=frame2)
    
    lab4 <- glabel("Step Size",cont=frame2)
    
    stepsizespin <- gedit(text="1", width = 20, cont=frame2)
    
    frame3 <- gframe("Plotting", cont=MainWinGroup, horizontal=F)
    
    linecheck <- gcheckbox("Linesplot", checked = TRUE, use.togglebutton=FALSE, handler = NULL, action = NULL,
                           container = frame3)
    
    linelegcheck <- gcheckbox("Line Plot Legend", checked = TRUE, use.togglebutton=FALSE, handler = NULL, action = NULL,
                              container = frame3)
    
    bpcheck <- gcheckbox("Label x axis by BP", checked = TRUE, use.togglebutton=FALSE, handler = NULL, action = NULL,
                         container = frame3)
    
    lab5 <- glabel("Frequency of x axis BP labelling", cont=frame3)
    bpspin <- gedit(text="800", width = 20, cont=frame3)
    
    mosaiccheck <- gcheckbox("Include Rainbow Bars", checked = TRUE, use.togglebutton=FALSE, handler = NULL, action = NULL,
                             container = frame3)
    
    lab6 <- glabel("Rainbow Bars Resolution", cont=frame3)
    mosaicspin <- gedit(text="200", width = 20, cont=frame3)
    
    combinecheck <- gcheckbox("Combine Plots?", checked = TRUE, use.togglebutton=FALSE, handler = NULL, action = NULL,
                              container = frame3)
    
    condensecheck <- gcheckbox("Condense Mosaics?", checked = TRUE, use.togglebutton=FALSE, handler = NULL, action = NULL,
                               container = frame3)
    
    lab7 <- glabel("Label font Sizes", cont=frame3)
    fontspin1 <- gedit(text="14", width = 20, cont=frame3)
    
    
    lab8 <- glabel("Legend font Sizes", cont=frame3)
    fontspin2 <- gedit(text="14", width = 20, cont=frame3)
    
    lab9 <- glabel("Specific Analysis", cont=frame3)
    analysistoplot <- gedit(text="NA", width=10, cont=frame3)
    
  } 
  
  .GlobalEnv$GUIconf <- conf
  
}