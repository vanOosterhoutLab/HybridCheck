---
layout: manualsection
title: Viewing and altering settings
manual: true
permalink: 02-manual.html
---

A word on HybRIDS various analysis settings
---------------------------------------------

This manual section deals with how to view and alter the HybRIDS analysis settings. But Before we do that we need to cover some aspects of the HybRIDS analysis and the settings it requires:

The HybRIDS analysis is split into several steps:

1. Triplet Generation
2. Sequence Similarity Analysis
3. Block Detection
4. Block Dating
5. Plotting

HybRIDS keeps the settings organised according to each of these steps.

---

Viewing analysis settings
-------------------------

**With the R console**

You can investigate a HybRIDS object by typing out its name which invokes it's show() method, this method prints to screen a summary of the object and all the information it contains:

```R
> MyAnalysis

HybRIDS object - Analysis of  3  aligned sequences.

DNA Alignment:
--------------
Full Sequence File Location:  /var/folders/kp/clkqvqn9739ffw2755zjwy74_skf_z/T//RtmproH4Bc/FullSequencec3bd78f6cb65 
Informative Sequence File Location:  /var/folders/kp/clkqvqn9739ffw2755zjwy74_skf_z/T//RtmproH4Bc/InformativeSequencec3bd12a391ee
Full Length:  398508 
Informative Length:  6577 
Sequence names:  Seq1 Seq2 Seq3 

Triplet Generation Parameters:
------------------------------
Triplet Generation Method: 1
Threshold for method number 2: 0.01

Sequence Similarity Analysis Parameters:
----------------------------------------
Sliding Window Size: 100
Sliding Window Step Size: 1

Block Detection Parameters: 
---------------------------
Manual Thresholds: 90
HybRIDS will attempt automatic detection of SS Thresholds for putative block searches.
HybRIDS will fall back on user specified manual thresholds, should the autodetection fail.

Block Dating Parameters:
------------------------
Assumed mutation rate: 1e-07
P-Value for acceptance of putative blocks: 0.005

No Triplets have been generated with the method makeTripletCombos yet.
```

You can see from the output, as we stated earlier - HybRIDS has settings associated with each of the 5 analysis steps listed previously.

So in short, to display the data and settings currently loaded into your HybRIDS object, as well as warning messages and other information, just type its name in the R console and hit enter.

You can also use a command `showParameters` on the HybRIDS object to get it to show the settings of a specific analysis step:

1. Specify "TripletGeneration" to see the settings for the triplet generation step: `showParameters("TripletGeneration")`
2. Specify "SSAnalysis" to see the settings for the Sequence Similarity Analysis step: `showParameters("SSAnalysis")`
3. Specify "BlockDetection" to see the settings for the Block Detection step: `showParameters("BlockDetection")`
4. Specify "BlockDating" to see the settings for the Block Dating step: `showParameters("BlockDating")`
5. Specify "Plotting" to see the settings for the plotting step: `showParameters("Plotting")`

As an example:

```R
> MyAnalysis$showParameters("BlockDetection")
Parameters for the detection of putative recombination blocks are:

A vector containing Manual Sequence Similarity Thresholds (%),

Whether or not you want HybRIDS to autodetect the sequence similarity thresholds,
                                                   
Whether you want HybRIDS to fall back and rely on the manual thresholds should the threshold autodetection fail.
                                                   
and finally a value by which the Standard Deviation of all sequence similarity is divded by during threshold detection (lower values = more conservative detection).

They are printed below.

$ManualThresholds
[1] 90

$AutoThresholds
[1] TRUE

$ManualFallback
[1] TRUE

$SDstringency
[1] 2

```

Many of HybRIDS settings have names that are self explnatory, but the manual page "HybRIDS Settings in Detail"
describes all of them in depth and the more commonly encountered/adjusted ones will be mentioned in the rest of the manual where applicable.

Note in the output of the last example of the `showParameters` method, a brief description of the settings is printed to console, before the values of the settings are printed.

**With the GUI**

With the GUI, a similar report is printed if you click the button "HybRIDS session summary button" which is located in the panel called "Data and Plotting".

---

Altering Analysis Settings
--------------------------

**With the R console**

Setting parameters in HybRIDS is done by using the method `setParameters` on a HybRIDS object. Two things need to be provided:

1. The name of the Analysis Step.
2. The names of the setting(s) you want to alter, assigned to their new value. 

To use an example, the line below specifies that HybRIDS should change the settings of the SSAnalysis (Sequence Similarity Analysis) step (aka step 2.), and then a setting with the value it should be changed to:

`MyAnalysis$setParameters("SSAnalysis", WindowSize = 50)`

This example sets the WindowSize variable of the SSAnalysis step to 50 base pairs.

This is how all settings in HybRIDS are set, note using this method you change multiple settings from the same analysis step at once:

`MyAnalysis$setParameters("BlockDating", MutationRate = 10e-9, DateAnyway = TRUE, PValue = 0.0005)`

This example changes three settings in the Block Dating analysis step.

**With the GUI**

With the GUI the manipulation of settings is a natural part of the usage of the interface:
When you click on one of the 5 analysis step buttons:

1. "Step 2). Generate the Triplet Sets for SSAnalysis".
2. "Step 3). Sequence Similarity Analysis".
3. "Step 4). Detect Putative Blocks".
4. "Step 5). Test and Date Putative Blocks".

A window will appear in which you can alter the settings before hitting the button to confirm and execute that step of the analysis.

---


