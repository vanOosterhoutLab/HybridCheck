<h1>Data and Code used in testing HybRIDS for publication.</h1>

This directory contains all data and code used in testing the HybRIDS package for the paper submitted to
Molecular Ecology Resources. It is a packrat based R project so all dependencies and requirements are stored in this project for maximum reproducibiliy of results.

This directory should be provided with the HybRIDS package in the extdata directory.

Here we list each file and describe what it contains / shows.

<h2>packrat</h2>
This folder contains the source of the R library and version of HybRIDS used to do this analysis.

<h2>Testing Power and Dating Folder:</h2>

<h3>Code Files</h3>

<dt><strong>Simulation_Functions.R</strong></dt>
<dd>Contains functions used to simulate sequence triplets containing recombination and process results. It is used by the files <code>Analyzing_Known_Events.R</code> and <code>Analyzing_Detected_Events</code>.</dd>
<dt><strong>Analyzing_Known_Events.R</strong></dt>
<dd>Simulates sequence triplets with a recombination breakpoint as described in the main text of the publication.
Then analyzes them with HybRIDS, informing HybRIDS of the locations of the recombinant blocks. This acts as a test of HybRIDS dating algorithm when all of the block is known. This code produces the datasets <code>KnownDatedAll.csv</code>, <code>KnownDatedSummary.csv</code>, and <code>SimulatedWithKnown.RData</code>.</dd>
  
<dt><strong>Analyzing_Detected_Events.R</strong></dt>
<dd>Simulates sequence triplets with a recombination breakpoint as described in the main text of the publication.
Then analyzes them with HybRIDS, first detecting blocks, and then dating the regions that are detected.
The code in this files produces the datasets <code>DetectionDatedAll.csv</code>, <code>DetectionDatedSummary.csv</code>, <code>detectionPercentagesFull.csv</code>, and <code>detectionPercentagesSummary.csv</code>, 
and <code>SimulatedWithDetection</code>.</dd>
</dl>
  
<h3>Data Files</h3>
<dl>
<dt><strong>SimulatedWithKnown.RData</strong></dt>
<dd>The raw blocks table output from simulated triplets tested and dated with HybRIDS when HybRIDS was given the locations of the blocks. This file was processed in <code>Analyzing_Known_Events.R</code> to produce <code>KnownDatedAll.csv</code> and
<code>KnownDatedSummary.csv</code>.</dd>

<dt><strong>SimulatedWithDetection.RData</strong></dt>
<dd>The raw blocks table output from simulated triplets tested and dated with HybRIDS when HybRIDS was used to detect
blocks. This file was processed in <code>Analyzing_Detected_Events.R</code> to produce <code>DetectionDatedAll.csv</code>,
<code>DetectionDatedSummary.csv</code>, <code>detectionPercentagesFull.csv</code>, and <code>detectionPercentagesSummary</code>.</dd>
  
<dt><strong>KnownDatedAll.csv</strong></dt>
<dd>Full processed dataset of the blocks detected and dated in simulations, when HybRIDS was given the locations of 
the recombination regions. The dataset shows how the estimated ages of blocks relate to the true ages of blocks,
and the amount of divergence between the parental sequences in the triplets.</dd>
  
<dt><strong>KnownDatedSummary.csv</strong></dt>
<dd>Summary processed dataset of the blocks detected and dated in simulations, when HybRIDS was given the locations of 
the recombination regions. The dataset shows how the estimated ages of blocks relate to the true ages of blocks,
and the amount of divergence between the parental sequences in the triplets. It is a summary of <code>KnownDatedAll.csv</code>.</dd>
  
<dt><strong>DetectionDatedAll.csv</strong></dt>
<dd>Full processed dataset of the blocks detected and dated in simulations, when HybRIDS detected recombination regions.
The dataset shows how the estimated ages of blocks relate to the true ages of blocks, and the amount of divergence
between the parental sequences in the triplets.</dd>
  
<dt><strong>DetectionDatedSummary.csv</strong></dt>
<dd>Summary processed dataset of the blocks detected and dated in simulations, when HybRIDS detected recombination
regions.
The dataset shows how the estimated ages of blocks relate to the true ages of blocks, and the amount of divergence
between the parental sequences in the triplets. It is a summary of <code>DetectionDatedAll.csv</code>.</dd>
   
<dt><strong>detectionPercentagesFull.csv</strong></dt>
<dd>Full processed dataset of the percentage of the recombination regions correctly detected by HybRIDS and how it relates to the amount of divergence between parental sequences, and the amount of divergence after the recombination
event in the triplets.</dd>
  
<dt><strong>detectionPercentagesSummary.csv</strong></dt>
<dd>Summary processed dataset of the percentage of the recombination regions correctly detected by HybRIDS and how it
relates to the amount of divergence between parental sequences, and the amount of divergence after the recombination
event in the triplets. It is a summary dataset of the dataset in <code>detectionPercentagesFull.csv</code>.</dd>
</dl>
  
<h2>Testing False Positive Rate Folder:</h2>

This folder contains code that was used to generate null scenarios - that is - triplets that don't show recombination.
These triplets were used to assess the false positive rate of HybRIDS and how it relates to the amount of divergence between the populations.

<h3>Code Files</h3>
<dl>
<dt><strong>Founding_Pop.py</strong></dt>
<dd>The simuPOP program that was used to generate the simulated sequences of populations that evolved in isolation.
It was run on a cluster and outputed CSV and GENEPOP files which were converted into fasta files with a Perl script.
The produced fasta files were the files used in HybRIDS.</dd>

<dt><strong>simuCSVtoFASTA2.pl</strong></dt>
<dd>This Perl script is used to generate FASTA formatted sequence files from the CSV files generated by the simuPOP
python script.</dd>
  
<dt><strong>False_Detection_Analysis.R</strong></dt>
<dd>Generates the files <code>fullFalsePos.csv</code> and <code>summaryFalsePos.csv</code>. It analyzes many triplets in which there is no
recombination, and calculates a false postive rate. The results allow one to see how the false positive rate is
affected by the amount of time the three sequences in the triplet have evolved independently of each other.</dd>
</dl>
    
<h3>Data Files</h3>
<dl>
<dt><strong>FASTA Sequences.zip</strong></dt>
<dd>This folder contains the FASTA formatted sequences output from <code>Founding_Pop.py</code> and converted into FASTA from 
CSV, by <code>simuCSVtoFASTA2.pl</code>.</dd>
  
<dt><strong>fullFalsePos.csv</strong></dt>
<dd>A data frame demonstrating how the false positive rate is affected by the amount of time the three sequences in the
triplet have evolved independently of one another.</dd>
  
<dt><strong>summaryFalsePos.csv</strong></dt>
<dd>A data frame that is a summary of <code>fullFalsePos.csv</code>. It also demonstrates how the (mean) false positive rate is
affected by the amount of time the three sequences in the triplet have evolved independently of one another.</dd>
</dl>
  

  



