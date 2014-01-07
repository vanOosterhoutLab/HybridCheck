---
layout: manualsection
title: Plotting and Data Export
manual: true
permalink: 05-manual.html
---

After all the analysis steps have been done, there are two more things that can be done - data export/tabulation, and plotting.

Data Export / Tabulation
------------------------

To collect the blocks data HybRIDS has produced and organize it into and R data frame, the tabulateDetectedBlocks() method is required.
Let's return to the MyAnalysis example and print out the blocks detected:

```R
table <- MyAnalysis$tabulateDetectedBlocks("Seq1:Seq2:Seq3")
```

A table of detected blocks for the only triplet present "Seq1:Seq2:Seq3" is assigned to the variable called 'table'. 
Look at such data tables by calling it's name on the command line:

```R
> table
   Sequence_Pair Sequence_Similarity_Threshold        Triplet First_BP_Position Last_BP_Position Approximate_Length_BP p=0.05_Age p=0.5_Age p=0.95_Age      P_Value
1      Seq1:Seq2                            89 Seq1:Seq2:Seq3            168661           170333                  1673       4900     16000      37600 4.444280e-05
2      Seq1:Seq2                            69 Seq1:Seq2:Seq3            170336           176546                  6211      22700     33300      46800 4.939782e-08
3      Seq1:Seq2                            69 Seq1:Seq2:Seq3            293405           296568                  3164      14900     27400      45600 8.383713e-06
4      Seq1:Seq2                            69 Seq1:Seq2:Seq3            366791           377186                 10396      20000     27600      37000 4.167285e-15
5      Seq1:Seq2                            69 Seq1:Seq2:Seq3            377388           378450                  1063       7700     25200      59200 4.376883e-03
6      Seq1:Seq3                            87 Seq1:Seq2:Seq3             51798            54807                  3010       4600     12200      25800 4.307145e-08
7      Seq1:Seq3                            49 Seq1:Seq2:Seq3             41449            51168                  9720      29300     38800      50100 4.922869e-08
8      Seq1:Seq3                            49 Seq1:Seq2:Seq3             54821            58301                  3481      15600     27800      45100 1.640235e-05
9      Seq1:Seq3                            49 Seq1:Seq2:Seq3             70217            91705                 21489      26100     32000      38800 3.743732e-21
10     Seq1:Seq3                            49 Seq1:Seq2:Seq3            239392           242595                  3204      10300     20900      37000 1.848867e-06
11     Seq1:Seq3                            49 Seq1:Seq2:Seq3            251167           271858                 20692      22800     28400      34900 9.343333e-24
12     Seq1:Seq3                            49 Seq1:Seq2:Seq3            383452           392106                  8655      24000     33200      44400 3.616045e-09
13     Seq2:Seq3                            81 Seq1:Seq2:Seq3            205714           219841                 14128       8300     12600      18100 7.453537e-32
14     Seq2:Seq3                            68 Seq1:Seq2:Seq3             26668            30534                  3867      12200     22500      37300 2.591180e-07
15     Seq2:Seq3                            68 Seq1:Seq2:Seq3            202629           205703                  3075       2700      8700      20500 1.777983e-09
```

Again as with most methods for HybRIDS objects, when you select a triplet(s) you can provide a vector of triplets e.g.
c("Seq1:Seq2:Seq3","Seq1:Seq2:Seq4....) or alternatively you can provide the shorthand version to choose all triplets
that contain the sequence pair you are interested in e.g. "Seq1:Seq2". In addition, this method has the optional
argument Neat which is TRUE by default. This clears up the table headers and omits some info the program uses and 
the user may be interested in but it's not needed, nor is it neat.

---

Plotting
--------

Plotting of sequence similarity values in HybRIDS is carried out by the plotSS() method.

Let's plot the first and only triplet from the MyAnalysis example from earlier:

```R
> MyAnalysis$plotSS("Seq1:Seq2:Seq3")
Selection Seq1:Seq2:Seq3 
```

This results in the graphic below.

![Can't Display Image](./img/exampleplot.jpeg)

The bottom graph is a line graph of sequence similarity between the three sequence pairs in the triplet i.e.
The sequence similarity between Seq1:Seq2, Seq1:Seq3, and Seq2:Seq3 for all the sliding window positions.
Areas where the values are high indicates a region of highly similar sequence that could be a recombination block.
The top graph shows the three sequences as coloured bars, Seq1 is red, Seq2 is green and Seq3 is blue.
In areas where the colours are the same and mixed indicates a high SS value and possible recombination region.
I.e. in a region where Seq1 and Seq2 are highly similar and there is a possible recombination event, the Red bar 
for Seq1 turns yellow as does the Green bar of Seq2 (red mixed with green = yellow in RGB colour mixing).
The MyAnalysis example sequences show quite a pronounced mosaic structure, with many "hybrid colour regions" 
highlighting possible recombinant events between the pairs of the sequence triplet.

---

**Plotting sequence similarity results for more than one HybRIDS Triplet**

Many triplets can contain a given pairing of two sequences, for example we know from the MyAnalysis2
triplets with pair "Seq1:Seq2" that there are several triplets containing Seq1 and Seq2.
It is useful therefore to compare the sequence similarity results between any pair of sequences,
of every triplet in which the pair occurs. We can do this with the same method as plotting a single triplet above,
but by specifying "Seq1:Seq2" instead of "Seq1:Seq2:seq3" to plot all sequence similarity results between Seq1 and Seq2
for all triplets in which they occur:

```R
> MyAnalysis2$plotSS("Seq1:Seq2")
Selection Seq1:Seq2
```

This takes a little bit of time to fetch all the data from all the triplets and transform it to make the data-frames needed to generate the plots. The result is a plot like the one below.

![Can't Display Image](./img/multiplotexample.jpeg)

We can see from this graph some interesting sequence regions where a highly dissimilar region between Seq1 and Seq2 is present in all triplets, and some highly similar regions occuring just downstream of that area for some triplets.

HybRIDS aims to provide many plotting parameters and options whilst providing reasonable and attractive defaults. There are many plotting parameters and so they are all described on the manual page describing all HybRIDS settings. The idea is to be able to quickly and easily change attributes of the plot without getting overly stuck into the R plotting environment and plotting libraries (of course the habitual R user knows exactly how to get the data out of the HybRIDS object and make their own custom plots).

When you are plotting, you can change the settings as you do for all other HybRIDS settings, but there is another way provided - but only for plotting:

You can change the plotting settings on the fly when you call the plotSS() method by feeding it extra parameters at the end of the function like so:

```R
#Let's plot only the bars:
> MyAnalysis$plotSS("Seq1:Seq2:Seq3", What="Bars")
Selection Seq1:Seq2:Seq3
``` 

You can also change the plotting parameters on the fly like above, but if you specify ReplaceParams as FALSE, then the change will only be temporary like so:

```R
> MyAnalysis$plotSS("Seq1:Seq2:Seq3", ReplaceParams=FALSE, What="Bars")
Selection Seq1:Seq2:Seq3
```

By setting ReplaceParams to false this allows you to experiment a few times with different plotting parameters without having to do the extra step of setting them with the setParameters() method all the time, and then setting them back if you don't like the result.

--- 








