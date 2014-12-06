---
layout: manualsection
title: Specifying sequence triplets
permalink: 05-manual.html
manual: true
published: true
status: publish
---
 
The first thing you need to do after creating the HybRIDS object for your alignment file, is inform the HybRIDS object which triplets for your sequences you would like to analyze for evidence of
recombination.
 
The following methods of generating triplets are provided, and are selected by changeing the `Method` parameter:
 
---
 
 **1.** (DEFAULT) - HybRIDS will analyze sequence triplets based on a set of specified groups (defined by the `Groups` parameter, such that triplets analyzed will be made up of sequences to try and find recombination events between the groups.
 
A group is any given subset of sequences from your input file. For example, the example anaysis has 10 DNA sequences. The first 3 sequences, might be from one population, the next three from a second population, and the remaining sequences from a third population. 
 
If you wanted to find evidence of recombination between these populations. You would want to analyze triplets made of e.g. The first sequence, the fourth sequence, and the seventh sequence, or e.g. the second sequence, the fourth sequence, and the ninth sequence. With this default method HybRIDS will figure out all the possible triplets that need to be analyzed that might find evidence of recombination between the populations or Groups. 
 
If this method of triplet generation is set, as it is by default, and no groups are provided, the HybRIDS object will analyze every possible combination of three of your provided sequences.
 
---
 
**2.** HybRIDS will analyze sequences triplet based on how similar the sequences that make the triplet are to one another. Every possible triplet of your provided sequences will be analyzed, providing that all the pairwise distances in a triplet are above a set threshold, and this threshold is defined by the `DistanceThreshold` parameter.
 
---
 
**3.** HybRIDS will analyze triplet based on sequence similarity as in method 2, however instead of using the `DistanceThreshold` parameter, a threshold is automatically chosen based on the distribution of pairwise distances for your provided alignment.
 
---
 
Now you know what this step does and how each of the settings for this step affect it, the example demonstrates how to appropriately edit the settings of this step.
 
In this example, there is an alignment of 10 sequences loaded in the HybRIDS object, and a researcer has collected the sequences from three populations/locations. Sequences one to three are from one location, sequences four to six are from the second location, and sequences 7 to 10 are from the final population or location.
 
First the sequences are loaded into a HybRIDS object as demonstrated in previous sections and the summary is printed to the console:
 

    library(HybRIDS)
    MyAnalysis <- HybRIDS$new("~/Dropbox/MySequences.fas")

    ## File to be read is expected to be FASTA format...
    ## Reading in sequence file...
    ## Looking for duplicates (sequences with p_distances of 0)...
    ## Done...
    ## Subsetting the informative segregating sites...
    ## Finished DNA input.
    ## Generating triplets to find recombination in sequences, between partitions.
    ## Initializing new triplets data.

    MyAnalysis

    ## HybRIDS object:
    ## 
    ## DNA Sequence Information:
    ## -------------------------
    ## An alignment of 10 sequences.
    ## 
    ## Full length of alignment: 400000
    ## Excluding non-informative sites: 33043
    ## 
    ## Sequence names:
    ## 1: Seq1
    ## 2: Seq2
    ## 3: Seq3
    ## 4: Seq4
    ## 5: Seq5
    ## 6: Seq6
    ## 7: Seq7
    ## 8: Seq8
    ## 9: Seq9
    ## 10: Seq10
    ## 
    ## Settings for Sequence Scan Combinations:
    ## ----------------------------------------
    ## Triplet Generation Method (Method): 1
    ## 
    ## Sequences are organized according to the following groups: 
    ## 
    ## 
    ## Distance Threshold for excluding triplets with
    ## 	too similar sequences (DistanceThreshold): 0.01
    ## 
    ## How many sequences from the same partition are
    ## 	allowed in a triplet (PartitionStrictness): 2
    ## 
    ## A total of 120 triplets will be compared.
    ## 
    ## A total of

    ## Error: object 'BlockDetectionParams' not found
 
We can see from the print-out that 120 triplets will be analyzed - this is because we have speified no groups but are using method 1. As a result, ever possible triplet that could be generated for these sequences will be analyzed.
 
Now the settings can be changed to add the groups between which we wish to find recombination, and edit the `PartitionStrictness` so as only one sequence from each group is allowed in a triplet.
 
In order to specify a group, you simply make a vector of the sequence names for each group, and then combine each one into a list:
 

    # Let's make our groups:
     
    FirstGroup <- c("Seq1", "Seq2", "Seq3")
    SecondGroup <- c("Seq4", "Seq5", "Seq6")
    ThirdGroup <- c("Seq7", "Seq8", "Seq9", "Seq10")
    myGroups <- list(FirstGroup, SecondGroup, ThirdGroup) # Put the groups in a list.
     
    # Let's edit the Triplet Generation settings
    MyAnalysis$setParameters("TripletGeneration", Groups = myGroups, PartitionStrictness = 1L)

    ## Generating triplets to find recombination in sequences, between partitions.
    ## Deleting all triplets data.
    ## Initializing new triplets data.

    # All done, now let's print the summary of our HybRIDS object again:
    MyAnalysis

    ## HybRIDS object:
    ## 
    ## DNA Sequence Information:
    ## -------------------------
    ## An alignment of 10 sequences.
    ## 
    ## Full length of alignment: 400000
    ## Excluding non-informative sites: 33043
    ## 
    ## Sequence names:
    ## 1: Seq1
    ## 2: Seq2
    ## 3: Seq3
    ## 4: Seq4
    ## 5: Seq5
    ## 6: Seq6
    ## 7: Seq7
    ## 8: Seq8
    ## 9: Seq9
    ## 10: Seq10
    ## 
    ## Settings for Sequence Scan Combinations:
    ## ----------------------------------------
    ## Triplet Generation Method (Method): 1
    ## 
    ## Sequences are organized according to the following groups: 
    ## c("Seq1", "Seq2", "Seq3"),
    ## c("Seq4", "Seq5", "Seq6"),
    ## c("Seq7", "Seq8", "Seq9", "Seq10")
    ## 
    ## Distance Threshold for excluding triplets with
    ## 	too similar sequences (DistanceThreshold): 0.01
    ## 
    ## How many sequences from the same partition are
    ## 	allowed in a triplet (PartitionStrictness): 1
    ## 
    ## A total of 36 triplets will be compared.
    ## 
    ## A total of

    ## Error: object 'BlockDetectionParams' not found
Note in the above example, that `1L` is used and not `1` when setting the `PartitionStrictness` as `1L` specifies the value is a whole number integer and not the floating point `1.0`.
 
We see after altering these settings, less triplets will be analyzed now, and each one will contain only one sequence from each group.
 
Now that you have instructed the HybRIDS object on which triplets will be analyzed, the recombination signal in each triplet can be scanned. This is explained in the next section.
---
 
 Step 2: Sequence Similarity Analysis
 ------------------------------------
 
 Once triplets are generated this step can be started by using the `analyzeSS()` method of the HybRIDS object.
 As a reminder the parameters of the step are:
 
 ```R
 > boop$showParameters("SSAnalysis")
Parameters for the Sliding Window Sequence Similarity analysis are:
Sliding Window Size,
Window Step Size, and the Sequence Triplet Combinations.
They are printed below.
 
$WindowSize
[1] 100
 
$StepSize
[1] 1
 
$TripletCombinations
$TripletCombinations[[1]]
[1] 1 2 3
```
 
Again these settings are described in the previous manual page, but by adjusting the WindowSize and Step size you adjust the size (in base pairs),
and the number of steps (in base pairs) that the sliding window slides along the Triplets during the analysis.
 
The SSAnalysis step is run very simply with the `analyzeSS()` method:
 
```R
> MyAnalysis$analyzeSS()
Only one triplet to analyze the sequence similarity of...
Preparing input DNA sequences...
Checking the sliding window parameters
Making all the window frames...
Analysing Now!
  |=============================================================================================================| 100%
``` 
 
Without any arguments, the `analyzeSS()` method will run the analysis for all triplets.
 
The analyzeSS method will analyse all the triplets by default, but you can also specify a triplet(s) to analyse:
This is useful if you have perhaps run the method for all your triplets (assuming you have more than one) and 
you've found few triplets in which you might try different window and step sizes to see if you can get different 
or better results. Here is a 10 Sequence example called "MyAnalysis2" to illustrate how this is done:
 
```
#Create the new HybRIDS object for a fasta file with 10 aligned sequences.
> MyAnalysis2 <- HybRIDS$new("~/Desktop/10Sequences.fasta")
 
Looking for duplicates...
 
#Make the triplet combinations or SSAnalysis step.
> MyAnalysis2$makeTripletCombos()
#Check out the MyAnalysis2 object before continuing.
> MyAnalysis2
HybRIDS object - Analysis of  10  aligned sequences.
 
DNA Alignment:
--------------
Full Sequence File Location:  /var/folders/kp/clkqvqn9739ffw2755zjwy74_skf_z/T//Rtmpqg9BMs/FullSequencec3d33691c805 
Informative Sequence File Location:  /var/folders/kp/clkqvqn9739ffw2755zjwy74_skf_z/T//Rtmpqg9BMs/InformativeSequencec3d31ce2dc40
Full Length:  400000 
Informative Length:  33043 
Sequence names:  Seq1 Seq2 Seq3 Seq4 Seq5 Seq6 Seq7 Seq8 Seq9 Seq10 
 
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
 
 120 Have been generated for analysis.
```
 
Now the example is set up, there's a 10 sequence (120 Triplet) HybRIDS object ready for sequence similarity analysis.
 
There are several methods in HybRIDS which can take a triplet selection as an option. 
To specify a triplet selection in HybRIDS, type out the sequence names, separated by semi-colons like so:
`"Seq1:Seq2:Seq3"`. This will select our first triplet in the MyAnalysis2 example. Note the order does not matter,
it could be `"Seq2:Seq1:Seq3"` or even `"Seq3:Seq2:Seq1"`.
 
So let's analyse only the Triplet `"Seq1:Seq2:Seq3"`. This is done like so:
 
```R
> MyAnalysis2$analyzeSS("Seq1:Seq2:Seq3")
Now analysing sequence similarity of triplet 1 2 3 
```
 
What if you want to do more than one selection? This is done by suppling a vector of several selections,
for example: `MyAnalysis2$analyzeSS(c("Seq1:Seq2:Seq3","Seq1:Seq2:Seq4","Seq1:Seq2:Seq5"))`.
 
Check the R documentation for use of the c() function - in short it is a concatenation function and you use 
it when you want to join several things into a vector - in this case we make a vector of our three triplet selections 
to give to HybRIDS.
 
There is an alternative shorthand to specifying Triplets in HybRIDS: say you wanted to analyze all the triplets that contained the two sequences called "Seq1" and "Seq2". You could do this:
 
```R
> MyAnalysis2$analyzeSS(c("Seq1:Seq2:Seq3","Seq1:Seq2:Seq4","Seq1:Seq2:Seq5","Seq1:Seq2:Seq6","Seq1:Seq2:Seq7","Seq1:Seq2:Seq8","Seq1:Seq2:Seq9","Seq1:Seq2:Seq10"))
Now analysing sequence similarity of triplet 1 2 3 
Now analysing sequence similarity of triplet 1 2 4 
Now analysing sequence similarity of triplet 1 2 5 
Now analysing sequence similarity of triplet 1 2 6 
Now analysing sequence similarity of triplet 1 2 7 
Now analysing sequence similarity of triplet 1 2 8 
Now analysing sequence similarity of triplet 1 2 9 
Now analysing sequence similarity of triplet 1 2 10 
```
 
But that gets awful to type out very quickly. In HybRIDS the shorthand for the above is this: `"Seq1:Seq2"` This amounts to telling HybRIDS The triplet must contain Seq1 and Seq2, and the third sequence is not specified so do it for all triplets with those two sequences.
Once again you could type this in any order: `"Seq1:Seq2"`, `"Seq2:Seq1"`. To demonstrate with an example:
 
```R
> MyAnalysis2$analyzeSS("Seq1:Seq2")
Now analysing sequence similarity of triplet 1 2 3 
Now analysing sequence similarity of triplet 1 2 4 
Now analysing sequence similarity of triplet 1 2 5 
Now analysing sequence similarity of triplet 1 2 6 
Now analysing sequence similarity of triplet 1 2 7 
Now analysing sequence similarity of triplet 1 2 8 
Now analysing sequence similarity of triplet 1 2 9 
Now analysing sequence similarity of triplet 1 2 10
```
 
In the example, HybRIDS does the analysis for all triplets of the 120 that contain both the first sequence and the 
second sequence (called "Seq1" and "Seq2" respectively).
 
---
 
 
Step 3: Block Detection Step
----------------------------
 
Identification of putative recombination blocks in triplets in HybRIDS is done using the findBlocks method, continuing with our MyAnalysis example:
 
```R
> test$findBlocks()
Only one triplet to find the potential blocks in...
Using the autodetect thresholds method...
Deciding on suitable thresholds...
Now beginning Block Search...
```
 
Again before using this method you can view and change the settings for that step.
 
As with the analyzeSS method, you can specify which triplets to find the blocks in:
 
```R
test$findBlocks("Seq1:Seq2:Seq3")
Now finding potential blocks in triplet 1 2 3 
Using the autodetect thresholds method...
Deciding on suitable thresholds...
Now beginning Block Search...
```
 
if you only use two sequence names - HybRIDS will perform the task but for all the triplets containing that pair of sequences. To demonstrate this with the 10 sequence example - MyAnalysis2:
 
```R
> MyAnalysis2$findBlocks("Seq1:Seq2")
Now finding potential blocks in triplet 1 2 3 
Using the autodetect thresholds method...
Deciding on suitable thresholds...
Falling back to manual thresholds...
Now beginning Block Search...
 
Now finding potential blocks in triplet 1 2 4 
Using the autodetect thresholds method...
Deciding on suitable thresholds...
Falling back to manual thresholds...
Falling back to manual thresholds...
Now beginning Block Search...
 
Now finding potential blocks in triplet 1 2 5 
Using the autodetect thresholds method...
Deciding on suitable thresholds...
Now beginning Block Search...
 
Now finding potential blocks in triplet 1 2 6 
Using the autodetect thresholds method...
Deciding on suitable thresholds...
Falling back to manual thresholds...
Falling back to manual thresholds...
Now beginning Block Search...
 
Now finding potential blocks in triplet 1 2 7 
Using the autodetect thresholds method...
Deciding on suitable thresholds...
Now beginning Block Search...
 
Now finding potential blocks in triplet 1 2 8 
Using the autodetect thresholds method...
Deciding on suitable thresholds...
Now beginning Block Search...
 
Now finding potential blocks in triplet 1 2 9 
Using the autodetect thresholds method...
Deciding on suitable thresholds...
Falling back to manual thresholds...
Now beginning Block Search...
 
Now finding potential blocks in triplet 1 2 10 
Using the autodetect thresholds method...
Deciding on suitable thresholds...
Now beginning Block Search...
```
 
The default options for this step will direct HybRIDS to try and auto-detect suitable sequence similarity 
thresholds for block detection, but will also fallback to the manual thresholds specified by 
the user should it not find any, this is default as it does not interrupt workflow should you 
be trying to run this for many triplets - all that will happen is HybRIDS falls back to the manual 
thresholds and then fails or succeeds to find blocks at those thresholds. 
 
---
 
Step 4: Block Dating (and Significance Testing)
-----------------------------------------------
 
Testing the significance and dating the blocks is as easy as it is to detect them. To do this use the dateBlocks method as demonstrated below:
 
```R
> MyAnalysis$dateBlocks()
Only one triplet to date blocks in...
Now dating blocks
```
 
Once again as with the analyzeSS and findBlocks methods, by default all triplets will be done, but you can provide a selection:
 
```R
> MyAnalysis$dateBlocks("Seq1:Seq2:Seq3")
Now dating blocks in triplet 1 2 3 
Now dating blocks
```
 
To demonstrate with the 10 sequence example - MyAnalysis2:
 
```R
> MyAnalysis2$dateBlocks("Seq1:Seq2")
Now dating blocks in triplet 1 2 3 
Now dating blocksNow dating blocks in triplet 1 2 4 
Now dating blocksNow dating blocks in triplet 1 2 5 
Now dating blocksNow dating blocks in triplet 1 2 6 
Now dating blocksNow dating blocks in triplet 1 2 7 
Now dating blocksNow dating blocks in triplet 1 2 8 
Now dating blocksNow dating blocks in triplet 1 2 9 
Now dating blocksNow dating blocks in triplet 1 2 10 
Now dating blocks
```
 
This will only test the significance of, and date significant blocks, for all triplets that contain both Seq1 and Seq2.
 
In ending this section, the main workflow and analysis process of HybRIDS has been covered.
The next step is really to decide how to view, assess and output your results. This is what will be covered in the next section.
