---
layout: manualsection
title: Loading HybRIDS and Entering Sequence Data
permalink: 02-manual.html
manual: true
published: true
status: publish
---
 
Starting HybRIDS and Loading Data
------------------------------------
 
To get started fire up the R condole and enter:
 

    library(HybRIDS)

    ## Loading required package: ggplot2
    ## Loading required package: methods
    ## Loading required package: grid
    ## Loading required package: gridExtra
    ## Loading required package: png
    ## Loading required package: ape

    ## 
    ## TACTTTGTACCTAAGTATGCATTACGTTACGTTAGTAGCTGGACCTAGTAAATCGGA     
    ## ,--.  ,--.         ,--.   ,------. ,--.,------.   ,---.
    ## |  '--'  |,--. ,--.|  |-. |  .--. '|  ||  .-.  \ '   .-'
    ## |  .--.  | \  '  / | .-. '|  '--'.'|  ||  |  \  :`.  `-.
    ## |  |  |  |  \   '  | `-' ||  |\  \ |  ||  '--'  /.-'    |
    ## `--'  `--'.-'  /    `---' `--' '--'`--'`-------' `-----'
    ##           `---'
    ## ATGAAACATGGATTCATACG - Version 1.0 - CGACCTGGATCATTTAGCCT
    ## 
    ## 
    ## Hybridisation, Recombination and Introgression Detection
    ##                    and Dating Package.
    ## 
    ## -----------------=========****=========------------------
    ## Cite: TBD
    ## Licence: GPL (Like R and most packages).
    ## http://ward9250.github.io/HybRIDS/
    ## -----------------=========****=========------------------
 
You should see the above result in the console:
 
The package is now successfully loaded.
 
---
 
Loading in sequence data and starting a new HybRIDS analysis
------------------------------------------------------------
 
To get started using HybRIDS to analyse a set of sequences, you first have to create a new HybRIDS object. Do this either with the R console or the Graphical Interface.
 
**With the R console:**
 
The HybRIDS object will contain all data, options, and will carry out the analysis. So everything you do in HybRIDS is done by interacting with this object.
By providing a file-path to a FASTA format sequence file when the new object is created the DNA data will also be loaded and prepared automatically.
 

    MyAnalysis <- HybRIDS$new("~/Desktop/MySequence.fasta")

    ## File to be read is expected to be FASTA format...
    ## Reading in sequence file...

    ## Warning: cannot open file '/local/yrq12edu/Desktop/MySequence.fasta': No
    ## such file or directory

    ## Error: cannot open the connection
With the above line I have created a new HybRIDS object and assigned it to a variable "MyAnalysis" - that is to say I made a new HybRIDS object that is called "MyAnalysis" and now it is in the R workspace ready for us to use by calling it's methods.
 
The HybRIDS package has been loaded, and data has been read in and is now contained in a HybRIDS object, which also contains all analysis settings and steps.
