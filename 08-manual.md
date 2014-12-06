---
layout: manualsection
title: Testing and dating detected blocks
permalink: 08-manual.html
manual: true
published: true
status: publish
---
 

 
Having searched for blocks in the desired sequence triplets, the `BlockDating` step is executed using the `dateBlocks` method.
 
This step assigns a significance value to blocks based on the size of a block, the number of mutations observed, and the binomial probability distribution. Afterwards, 95%CI regions for the divergence time of the blocks are calculated, assuming a mutational model and a strict molecular clock.
 
The parameters for this step are:
 
MutationRate
  : The substitution rate assumed when estimating the ages of blocks. By default it is 10e-9.
 
PValue
  : The critical alpha value when testing the significance values of blocks. By default is is 0.05.
  
BonfCorrection
  : Logical (`TRUE` / `FALSE`) when `TRUE` (default) the `PValue` will be adjusted with a Bonferroni correction.
  
DateAnyway
  : Logical (`TRUE` / `FALSE`) when `TRUE` (`FALSE` by default), all detected blocks will be kept and dated.
  : When `FALSE`, blocks that do not pass the critical alpha when testing for significance are not kept in the results
  : and are not dated.
  
MutationCorrection
  : Set to one of "raw", "TS", "TV", "JC69", "K80", "F81", "K81", "F84", "BH87", "T92", "TN93", "GG95".
  : Base frequencies are calculated from the recombinant block if they are required by the model.
 
 
As was the case with the `analyzeSS` and `findBlocks` methods, `dateBlocks` needs a first argument that dictates which triplets of the set defined by the `TripletGeneration` parameters will have their blocks tested and dated. By default this is `"NOT.DATED"`, which means all triplets with blocks detected, that have not been tested or dated, will have their detected blocks tested and dated witht the current settings for the step. This can be set to `"ALL"`, or, as previously, can be set as a list of vectors that each contain three sequence names.
 
The following example finds blocks in the two triplets that were scanned for recombination signal and had blocks identified in the examples of the previous two sections:

    tripletsToDate <- list(c("Seq1", "Seq4", "Seq7"), c("Seq1", "Seq4", "Seq8"))
    MyAnalysis$dateBlocks(tripletsToDate)

    ## Now dating blocks
    ## Now dating blocks
 
The blocks can then be tabulated again, the columns of the table that were previously `NA` are now filled with values:
 

    MyAnalysis$tabulateDetectedBlocks(tripletsToDate)

    ## Tabulating blocks for the triplet Seq1:Seq4:Seq7
    ## Tabulating blocks for the triplet Seq1:Seq4:Seq8

    ##          Triplet Sequence_Pair Sequence_Similarity_Threshold
    ## 1 Seq1:Seq4:Seq7     Seq1:Seq4                            65
    ## 2 Seq1:Seq4:Seq7     Seq1:Seq7                            85
    ## 3 Seq1:Seq4:Seq7     Seq4:Seq7                            68
    ## 4 Seq1:Seq4:Seq8     Seq1:Seq4                            72
    ## 5 Seq1:Seq4:Seq8     Seq1:Seq8                            84
    ## 6 Seq1:Seq4:Seq8     Seq4:Seq8                            66
    ##   First_BP_Position Last_BP_Position Approximate_Length_BP Number_of_SNPs
    ## 1            287454           295896                  8443             56
    ## 2            283984           287383                  3400             21
    ## 3            269143           280337                 11195             71
    ## 4            287463           295890                  8428             55
    ## 5            283980           287386                  3407             22
    ## 6            269088           280341                 11254             75
    ##   Corrected_Number_of_SNPs p=0.05_Age p=0.5_Age p=0.95_Age   P_Value
    ## 1                       56      41443     33542      26885 3.877e-63
    ## 2                       21      44389     31743      21901 5.136e-39
    ## 3                       71      38663     31994      26187 1.491e-65
    ## 4                       55      41014     33011      26335 1.120e-63
    ## 5                       22      46162     33147      23016 1.526e-38
    ## 6                       75      40223     33579      27611 1.935e-63
    ##   P_Thresh
    ## 1    0.005
    ## 2    0.005
    ## 3    0.005
    ## 4    0.005
    ## 5    0.005
    ## 6    0.005
 
  
  
 
 
 
 
 
