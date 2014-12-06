---
layout: manualsection
title: Scanning triplets for recombination signal
permalink: 06-manual.html
manual: true
published: true
status: publish
---
 

 
This step can be started by using the `analyzeSS()` method of the HybRIDS object.
As a reminder the parameters of the step are:
 
WindowSize
  : The size in bases of the Sliding Window. The default size is 100.
 
StepSize
  : The number of bases that the sliding window moves across sequences. The default size is 1. 
 
The `analyzeSS` method also accepts an argument which tells it which triplets
(of the set you have defined in the `TripletGeneration` step) are to be scanned with the current
`SSAnalysis` settings.
 
This argument can be set as `"NOT.SCANNED"` (default) in which case, every triplet of the set defined by the 
`TripletGeneration` parameters will be scanned for recombination signal with the current `SSAnalysis`
settings, if it has not previously already been scanned.
This argument can also be set to `"ALL"`, in which case all the triplets defined by the `TripletGeneration`
settings will be scanned for recombination signal.
 
If the argument is not one of these two arguments, it can be a list of vectors that each contain three sequence names.
 
For example just two triplets could be scanned:
 

    tripletsToScan <- list(c("Seq1", "Seq4", "Seq7"), c("Seq1", "Seq4", "Seq8")) 
     
    MyAnalysis$analyzeSS(tripletsToScan)

    ## Scanning sequence similarity for triplet Seq1, Seq4, Seq7
    ## Checking the sliding window parameters...
    ## Making all the window frames...
    ## Scanning Now!
    ## Scanning sequence similarity for triplet Seq1, Seq4, Seq8
    ## Checking the sliding window parameters...
    ## Making all the window frames...
    ## Scanning Now!
    ## Finished Sequence Similarity Analysis.
 
 
 
 
 
 
 
 
 
