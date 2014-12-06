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
 
**1.**
	: (DEFAULT) - HybRIDS will analyze sequence triplets based on a set of specified groups (defined by the `Groups` parameter, such that triplets analyzed will be made up of sequences to try and find recombination events between the groups.
 
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
    ## Settings for sliding window scan of recombination signal:
    ## ---------------------------------------------------------
    ## Size of the sliding window in base pairs (WindowSize): 100
    ## 
    ## Number of base pairs to move the sliding window on
    ## 	each iteration of the scan (StepSize): 1
    ## 
    ## A total of
 
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
    ## Settings for sliding window scan of recombination signal:
    ## ---------------------------------------------------------
    ## Size of the sliding window in base pairs (WindowSize): 100
    ## 
    ## Number of base pairs to move the sliding window on
    ## 	each iteration of the scan (StepSize): 1
    ## 
    ## A total of
Note in the above example, that `1L` is used and not `1` when setting the `PartitionStrictness` as `1L` specifies the value is a whole number integer and not the floating point `1.0`.
 
We see after altering these settings, less triplets will be analyzed now, and each one will contain only one sequence from each group.
 
Now that you have instructed the HybRIDS object on which triplets will be analyzed, the recombination signal in each triplet can be scanned. This is explained in the next section.
 
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
