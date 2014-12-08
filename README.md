# Data and Code used in testing HybRIDS for publication.

This directory contains all data and code used in testing the HybRIDS package for the paper submitted to
Molecular Ecology Resources. It is a packrat based R project so all dependencies and requirements are stored in this project for maximum reproducibiliy of results.

This directory should be provided with the HybRIDS package in the extdata directory.

Here we list each file and describe what it contains / shows.

## Code:

**Simulation_Functions.R**
  : Contains functions used to simulate sequence triplets containing recombination and process results. It is used by
  : the files `Analyzing_Known_Events.R` and `Analyzing_Detected_Events`.
  
**Analyzing_Known_Events.R**
  : Simulates sequence triplets with a recombination breakpoint as described in the main text of the publication.
  : Then analyzes them with HybRIDS, informing HybRIDS of the locations of the recombinant blocks. This acts
  : as a test of HybRIDS dating algorithm when all of the block is known. This code produces the datasets
  : `KnownDatedAll.csv`, `KnownDatedSummary.csv`, and `SimulatedWithKnown.RData`.
  
**Analyzing_Detected_Events.R**
  : Simulates sequence triplets with a recombination breakpoint as described in the main text of the publication.
  : Then analyzes them with HybRIDS, first detecting blocks, and then dating the regions that are detected.
  : The code in this files produces the datasets `DetectionDatedAll.csv`, `DetectionDatedSummary.csv`,
  : `detectionPercentagesFull.csv`, and `detectionPercentagesSummary.csv`, and `SimulatedWithDetection`.
  
## Data Files

**SimulatedWithKnown.RData**
  : The raw blocks table output from simulated triplets tested and dated with HybRIDS when HybRIDS was given the
  : locations of the blocks. This file was processed in `Analyzing_Known_Events.R` to produce `KnownDatedAll.csv` and
  : `KnownDatedSummary.csv`.

**SimulatedWithDetection.RData**
  : The raw blocks table output from simulated triplets tested and dated with HybRIDS when HybRIDS was used to detect
  : blocks. This file was processed in `Analyzing_Detected_Events.R` to produce `DetectionDatedAll.csv`,
  : `DetectionDatedSummary.csv`, `detectionPercentagesFull.csv`, and `detectionPercentagesSummary`.
  
**KnownDatedAll.csv**
  : Full processed dataset of the blocks detected and dated in simulations, when HybRIDS was given the locations of 
  : the recombination regions. The dataset shows how the estimated ages of blocks relate to the true ages of blocks,
  : and the amount of divergence between the parental sequences in the triplets.
  
**KnownDatedSummary.csv**
  : Summary processed dataset of the blocks detected and dated in simulations, when HybRIDS was given the locations of 
  : the recombination regions. The dataset shows how the estimated ages of blocks relate to the true ages of blocks,
  : and the amount of divergence between the parental sequences in the triplets. It is a summary of `KnownDatedAll.csv`.
  
**DetectionDatedAll.csv**
  : Full processed dataset of the blocks detected and dated in simulations, when HybRIDS detected recombination regions.
  : The dataset shows how the estimated ages of blocks relate to the true ages of blocks, and the amount of divergence
  : between the parental sequences in the triplets.
  
**DetectionDatedSummary.csv**
  : Summary processed dataset of the blocks detected and dated in simulations, when HybRIDS detected recombination
  : regions.
  : The dataset shows how the estimated ages of blocks relate to the true ages of blocks, and the amount of divergence
  : between the parental sequences in the triplets. It is a summary of `DetectionDatedAll.csv`.
   
**detectionPercentagesFull.csv**
  : Full processed dataset of the percentage of the recombination regions correctly detected by HybRIDS and how it
  : relates to the amount of divergence between parental sequences, and the amount of divergence after the recombination
  : event in the triplets. 
  
**detectionPercentagesSummary.csv**
  : Summary processed dataset of the percentage of the recombination regions correctly detected by HybRIDS and how it
  : relates to the amount of divergence between parental sequences, and the amount of divergence after the recombination
  : event in the triplets. It is a summary dataset of the dataset in `detectionPercentagesFull.csv`.
  
## In the simuPOP simulations folder

This folder contains code that was used to generate null scenarios - that is - triplets that don't show recombination.
These triplets were used to assess the false positive rate of HybRIDS and how it relates to the amount of divergence between the populations.

**Founding_Pop.py**
  : The simuPOP program that was used to generate the simulated sequences of populations that evolved in isolation.
  : It was run on a cluster and outputed CSV and GENEPOP files which were converted into fasta files with a Perl script.
  : The produced fasta files were the files used in HybRIDS.
  
**FASTA Sequences**
  : This folder contains the CSV, and, FASTA formatted sequences output from `Founding_Pop.py`.
  
**simuCSVtoFASTA2.pl**
  : This Perl script is used to generate FASTA formatted sequence files from the CSV files generated by the simuPOP
  : python script.
  
**False_Detection_Analysis.R**
  : Generates the files `fullFalsePos.csv` and `summaryFalsePos.csv`. It analyzes many triplets in which there is no
  : recombination, and calculates a false postive rate. The results allow one to see how the false positive rate is
  : affected by the amount of time the three sequences in the triplet have evolved independently of each other.
  
**fullFalsePos.csv**
  : A data frame demonstrating how the false positive rate is affected by the amount of time the three sequences in the
  : triplet have evolved independently of one another.
  
**summaryFalsePos.csv**
  : A data frame that is a summary of `fullFalsePos.csv`. It also demonstrates how the (mean) false positive rate is
  : affected by the amount of time the three sequences in the triplet have evolved independently of one another.
  

  



