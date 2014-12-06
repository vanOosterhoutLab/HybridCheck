---
layout: manualsection
title: Detection of recombination blocks in signal
permalink: 07-manual.html
manual: true
published: true
status: publish
---
 

 
Having scanned for recombination signal for the desired sequence triplets, the `BlockDetection` step is executed using the `findBlocks` method.
 
This step identifies regions of recombination in a given sequence triplet by examining the data from the sliding window scans in the previous `SSAnalysis` step. If there are any regions in which two sequences meet a certain sequence similarity threshold, they are treated as potentially recombinant.
 
The parameters for this step are:
 
ManualThresholds
  : A vector of numbers between 0 and 100.
  : Each is a percenage sequence similarity threshold.
  : These thresholds are set by the user, and can be used to find recombination,
  : or used as a fallback if HybRIDS fails to automatically decide on some thresholds.
  
AutoThresholds
  : A logical argument (`TRUE`, or `FALSE`), if `TRUE` (default) then during execution of the `findBlocks` method,
  : thresholds will be decided on automatically. If `FALSE` then the `ManualThresholds` will be used.
  
ManualFallback
  : A logical argument (`TRUE`, or `FALSE`), if `TRUE` (default), then if automatic thresholds cannot be found, the
  : `ManualThresholds` will be used instead.
  
SDstringency
  : A numeric value, the higher it is, the closer automatically detected thresholds are allowed to be to the mean
  : sequence similarity across the sequences. It is usually fine to leave this as 2.
  
Like the `analyzeSS` method in the previous section, the `findBlocks` must be given an argument that dictates which triplets of the set defined by the `TripletGeneration` parameters will have blocks found from their scan data.
By default this is set to `"NOT.SEARCHED"` which will find blocks for triplets which have not already had blocks identified from their scan data. This argument can also be set to `"ALL"` or a list of vectors that each contain three sequence names.
 
The below example finds blocks in the two triplets that were scanned for recombination signal in the previous section:

    tripletsToSearch <- list(c("Seq1", "Seq4", "Seq7"), c("Seq1", "Seq4", "Seq8"))
    MyAnalysis$findBlocks(tripletsToSearch)

    ## Using the autodetect thresholds method...
    ## Deciding on suitable thresholds...
    ## Now beginning Block Search...
    ## Using the autodetect thresholds method...
    ## Deciding on suitable thresholds...
    ## Now beginning Block Search...
    ## Finished finding potential blocks for all triplet selections.
 
Note that, as with the `analyzeSS` method, you cannot specify a triplet which was not specified by the `TripletGeneration` settings:
 

    tripletsToSearch <- list(c("Seq1", "Seq2", "Seq3"))
    MyAnalysis$findBlocks(tripletsToSearch)

    ## Finished finding potential blocks for all triplet selections.
 
You also cannot search for blocks for triplets which have not been analyzed by `analyzeSS`:

    tripletsToSearch <- list(c("Seq1", "Seq4", "Seq9"))
    MyAnalysis$findBlocks(tripletsToSearch)

    ## Finished finding potential blocks for all triplet selections.

    # A triplet must be scanned for recombination signal first
    MyAnalysis$analyzeSS(tripletsToSearch)

    ## Finished Sequence Similarity Analysis.

    MyAnalysis$findBlocks(tripletsToSearch)

    ## Finished finding potential blocks for all triplet selections.
 
If you want to get a Data Frame containing a table of the detected blocks in some triplets, use the `tabulateDetectedBlocks` method.
 
This method needs to be given a set of sequence names to set the triplets that will be tabulated:
 

    tripletsToTabulate <- list(c("Seq1", "Seq4", "Seq7"), c("Seq1", "Seq4", "Seq7"))
    MyAnalysis$tabulateDetectedBlocks(tripletsToTabulate)

    ## Tabulating blocks for the triplet Seq1:Seq4:Seq7

    ##          Triplet Sequence_Pair Sequence_Similarity_Threshold
    ## 1 Seq1:Seq4:Seq7     Seq1:Seq4                            65
    ## 2 Seq1:Seq4:Seq7     Seq1:Seq7                            85
    ## 3 Seq1:Seq4:Seq7     Seq1:Seq7                            68
    ## 4 Seq1:Seq4:Seq7     Seq1:Seq7                            68
    ## 5 Seq1:Seq4:Seq7     Seq4:Seq7                            68
    ##   First_BP_Position Last_BP_Position Approximate_Length_BP Number_of_SNPs
    ## 1            287454           295896                  8443             NA
    ## 2            283984           287383                  3400             NA
    ## 3            283957           283980                    24             NA
    ## 4            287386           287407                    22             NA
    ## 5            269143           280337                 11195             NA
    ##   Corrected_Number_of_SNPs p=0.05_Age p=0.5_Age p=0.95_Age P_Value
    ## 1                       NA         NA        NA         NA      NA
    ## 2                       NA         NA        NA         NA      NA
    ## 3                       NA         NA        NA         NA      NA
    ## 4                       NA         NA        NA         NA      NA
    ## 5                       NA         NA        NA         NA      NA
    ##   P_Thresh
    ## 1       NA
    ## 2       NA
    ## 3       NA
    ## 4       NA
    ## 5       NA
 
 
 
 
