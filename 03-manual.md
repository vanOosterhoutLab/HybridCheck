---
layout: manualsection
title: Viewing and editing analysis settings
permalink: 03-manual.html
manual: true
published: true
status: publish
---
 

 
 
Viewing HybRIDS current settings
--------------------------------
 
To view the settings and state of your current HybRIDS object, you just need to enter the name of the variable you used to assign the HybRIDS object, in the case of this example, the HybRIDS object created was assigned to the name *MyAnalysis*
 

    MyAnalysis

    ## HybRIDS object:
    ## 
    ## DNA Sequence Information:
    ## -------------------------
    ## An alignment of 3 sequences.
    ## 
    ## Full length of alignment: 398508
    ## Excluding non-informative sites: 6496
    ## 
    ## Sequence names:
    ## 1: Seq1
    ## 2: Seq2
    ## 3: Seq3
    ## 
    ## Settings for Sequence Scan Combinations:
    ## ----------------------------------------
    ## Triplet Generation Method (Method): 1
    ## 
    ## Distance Threshold for eliminating triplets with
    ## 	too similar sequences (DistanceThreshold): 0.01
    ## 
    ## How many sequences from the same partition are
    ## 	allowed in a triplet (PartitionStrictness): 2
    ## 
    ## A total of 1 triplets will be compared.
    ## 
    ## A total of

    ## Error: object 'BlockDetectionParams' not found
You'll see a text summary printed like that above which shows the settings for each analysis step.
 
A full description of each setting is given in the next manual section.
