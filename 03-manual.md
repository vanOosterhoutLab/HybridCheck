---
layout: manualsection
title: Viewing and editing analysis settings
permalink: 03-manual.html
manual: true
published: true
status: publish
---
 

 
When scripting with a HybRIDS object, you execute several methods that make the HybRIDS object perform several steps of the analysis. Each analysis step is performed according to its settings. Before we walk through executing the analysis steps in detail, it's important to know how to view and alter the settings of the steps.
 
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
 
A full description of each setting is given in the upcoming manual sections dedicated to each step.
 
Changeing the analysis settings
-------------------------------
The settings of each analysis step is altered in a consistent way, using the same method of the HybRIDS object.
 
To change the settings for any given step of the HybRIDS analysis, you use the `setParameters` method.
 
This method accepts first a word that dictates which analysis step it is you want to change:
 
1. `TripletGeneration` - Settings for Sequence Scan Combinations.
2. `SSAnalysis` - Settings for the sequence similarity scan step.
3. `BlockDetection` - Settings for detecting recombinant blocks from sequence similarity scan data generated from the sequence similarity scan step.
4. `BlockDating` - Settings for calculating the significance values and divergence time estimates, for recombinant blocks detected in the `BlockDetection` step.
 
In addition the method is then provided a series of name, value pairs where a name is one of the settings, and the value is the value to change it to.
 
To use an example, we know from the text summary printed for the HybRIDS object above, that the `TripletGeneration` step has several settings: `Method`, `DistanceThreshold`, and `PartitionStrictness`.
 
(Don't worry about what these settings do right now, they are described more fully in the next manual sections - although the nice text summary gives you a clue e.g. `DistanceThreshold` is a *Distance Threshold for eliminating triplets with too similar sequences*)
 
Let's say we wanted to increace the `DistanceThreshold` parameters from 0.01 to 0.1. We can do that by calling the `setParameters` method:

    MyAnalysis$setParameters("TripletGeneration", DistanceThreshold = 0.1)

    ## Warning: There is only one comparrison possible - presumably only 3
    ## sequences are present.

    ## Deleting all triplets data.
    ## Initializing new triplets data.
 
Multiple settings could be edited too, for example:
```R
MyAnalysis$setParameters("TripletGeneration", DistanceThreshold = 0.1, PartitionStrictness = 1)
```
 
To summarize: Each step of a HybRIDS analysis has settings associated with them, but the settings of each step are all easily changed by calling the `setParameters` method of the HybRIDS object.
 
In the next sections, each step of the HybRIDS analysis, and how to execute them, will be described.
