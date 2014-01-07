---
layout: manualsection
title: Descriptions of all HybRIDS settings
manual: true
permalink: 03-manual.html
---

Triplet Generation Step Settings
--------------------------------

* **Method:**
This variable sets the method by which HybRIDS will generate the sequence triplets that will be analyzed. 
There are three methods that can be set 1, 2 or three, the methods are as follows:
	1. Doing the triplet generation step with this first method means that for the remainder of the analysis, every possible
	combination of three sequences (called a triplet) will be prepped for by HybRIDS - internal variables will be set accordingly 
	and placeholders for the appropriate number of triplets will be prepared to contain generated data.
	2. Selecting this method will use distance information is used to only prepare for analysing a set of triplets which do not contain sequence pairs
	which would likeley make for a poor analysis because of the insufficient genetic distance between them and inadequate number of informative base positions.
	The user can specify a threshold for this method by adjusting the value of the "SortThreshold" parameter of this analysis step.
	3. Selecting this method will also only prepare for analysing a set of triplets which do not contain sequence pairs
	which would likeley make for a poor analysis because of the insufficient genetic distance between them and inadequate number of informative base positions.
	However instead of a user defined variable, HybRIDS auto-decides on a threshold based on the distribution of distances between all sequence pairs, eliminating 
	those lower and apart from the rest of the distribution i.e. low outliers.

<br>

* **SortThreshold:** 
This variable sets the threshold used for the second method of triplet generation described above. 

---
	
Sequence Similarity Analysis Step Settings
------------------------------------------

* **WindowSize:**
This variable sets the size of the sliding window (in bp) that moves across sequence triplets, measuring the sequence similarities for that segment. 
Lower values mean smaller windows and a larger value means bigger windows. Larger windows mean more averaging of the data and so removal of noisy fluctuations but at the same time possible elimination of genuine recombination signal.
Smaler windows mean more noise and possibly more false positive but my be beneficial if you are looking for smaller regions or at smaller sequences.

* **StepSize:**
This variable sets the step size of the sliding window (in bp), this denotes the number of bp by which the sliding window moves along. A step size of one means many windows, but finer scale measurements of sequence similarity changes across sequence triplets.
Whereas larger step sizes mean bigger window jumps and so less windows in total and so you get a coarser view of sequence similarity changes across a sequence triplet, but you save on memory and time, which may be critical for extremely large datasets. 

* **TripletCombinations:**
This is intended mostly as an internal variable and so should be left alone.

---

Block Detection Step Settings
-----------------------------

* **ManualThresholds:**
This is a vector (in the R sense) of numbers that form Sequence Similarity Threshold percentage values for block detection. 
Any region of the sequence triplet in which two sequences are so similar that one of these are exceeded, is flagged as a potential recombinant block between the two sequences and data about the region is saved. You should not specify numbers out side the range 0 (i.e. 0% = totally different sequences) - 100 (i.e. 100% = totally identical sequences).

* **AutoThresholds:**
This is a boolean value set to either true for false. If true, HybRIDS will auto decide the thresholds to use for block detection by analysing the distribution
of sequence similarity values for the entire length of the sequence triplet, and identify multiple modes/peaks in the distribution and decide if those values stand apart from the noise of the distribution, if so it is treated as a threshold for block detection.

* **ManualFallback**
This is a boolean value set to either true or false. If true, after attmpting to auto-decide the thresholds for putative recombinant block detection, HybRIDS will - if it fails to auto-decide on some thresholds - fall back and use the thresholds specified by the **ManualThresholds** setting.

* **SDstringency**
This is a value which lessens or tightens the stringency of the auto-matic block detection process. To pass the putative block detection, the sequence similarity between two sequences must be higher than the mean sequence similarity level + the standard deviation divided by this value, by default it is two, which is equal to half the standard deviation.

---

Block Dating Settings
---------------------

* **MutationRate:**
This is the mutation rate that is to be supplied to the date estimation part of the algorithm, the closer those value is to the true mutation rate of the sequences under study the better.

* **PValue:**
This is the starting P-Value threshold to be provided to the putative block testing part of the algorithm: For each block putatively detected in the previous "Block Detection Step",
a significance test is applied which results in a P-Value which is the probability the region is result of mutation alone. A P-Value of 0.005 means a 0.5% change the region is result of mutation alone 
and so it is highly likely the result of recombination or similar process.
This value is adjusted by bonferroni correction and blocks with P-Values that are not smaller than the adjusted threshold are not dated or included in the final result output.

* **BonfCorrection:**
This is a boolean value that if set to false, turns off the bonferroni correction of the PValue parameter during significance testing.

* **DateAnyway:**
This is a boolean value that can be set to true or false, if set to true, even if a putative block region does not pass the significance test because it's P-value is higher than the threshold,
it will still get dated and included with the results.

---

Plotting Settings
-----------------

* **What:**
This value can either be "Bars", "Lines" or an R vector of both. This determines which plots HybRIDS will plot. HybRIDS produces two class of plot of a given triplet's SS similarity data (often on one canvas), a lines plot, and a coloured heatmap-like plot.

* **PlotTitle:**
Boolean, set wether to include a title with the plot.

* **CombinedTitle:**
Boolean, set whether to only have one overall title for a Lines and Bars plot on the same canvas. 

* **TitleSize:**
Integer value, sets the font size of the text in the plot title. 

* **TitleFace:**
A character value "bold" by default, can be changed to a face compatible with ggplot2, see documentation of ggplot2 for more details.

* **TitleColour:** 
A string denoting a colour "black" by default. Can be changed to any colour notation compatible with ggplot2, see the documentation or examples of ggplot2 for more details.

* **XLabels:**
Boolean, set whether to include the value labels of the x-axis in the plots. 

* **YLabels:**
Boolean, set whether to include the value labels of the y-axis in the plots.

* **XTitle:**
Boolean, set whether to include the title of the x-axis in the plots.

* **XTitleFontSize:**
Integer value, sets the font size of the x-axis title in plots.

* **XTitleColour:**
A character string denoting a colour, "black" by default. Sets the colour of the x-axis title font. Can be changed to any colour notation compatible with ggplot2, see the documentation or examples of ggplot2 for more details.

* **XLabelSize:**
Integer value, sets the font size of the x-axis labels.

* **XLabelColour:**
A character string denoting a colour, "black" by default. Sets the colour of the x-axis value labels. Can be changed to any colour notation compatible with ggplot2, see the documentation or examples of ggplot2 for more details.

* **YTitle:**
Boolean, set whether to include the title of the y-axis in the plots.

* **YTitleFontSize:**
Integer value, sets the font size of the y-axis title in plots.

* **YTitleColour:**
A character string denoting a colour, "black" by default. Sets the colour of the y-axis title font. Can be changed to any colour notation compatible with ggplot2, see the documentation or examples of ggplot2 for more details.

* **YLabelSize:**
Integer value, sets the font size of the y-axis labels.

* **YLabelColour:**
A character string denoting a colour, "black" by default. Sets the colour of the x-axis value labels. Can be changed to any colour notation compatible with ggplot2, see the documentation or examples of ggplot2 for more details.

* **Legends:**
Boolean value, set's whether to include the legends of plots. TRUE by default.

* **LegendFontSize:**
An integer value, sets the legend font size.

* **MosaicScale:**
An integer value, 500 by default. This affects how the Bars plot is drawn. With a higher value, a greater resolution is achieved but more and finer sequence similarity data is required. This is achieved with sequences with high numbers of informative sites and by adjusting the WindowSize and StepSize settings for the SSAnalysis step.
HybRIDS will print a message to the console if problems arise during plot drawing due to the MosaicScale settings.

---




