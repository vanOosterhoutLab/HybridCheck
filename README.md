HybRIDS
=======

**An R package for simple visualisation of recombinant regions and mosaic genome structures. Including block detection and simple dating estimation.**

# Introduction

HybRIDS is an R package that began life as some code to allow colleagues to attractively visualize the mosaic structure of
a sequence - what was wanted turned out to be more complex than usual plotting in R the casual user may usually do and required use of several libraries and a non-trivial amount of data-munging and code.

Therefore we decided to provide this package for other researchers so they can quickly visualize the mosaic structure of sequences in a visually appealing way - thanks to R's wonderful graphical system - 
along with the ability to detect potential blocks based on the analysis method, and also to date those putative blocks. We then also coded a graphical user interface for the package to assist the occasional R user.  
HybRIDS then can cater to the individual who is comfortable with code and who wants to script many analyses steps together and want's to include some HybRIDS visualisation and dating, and also to the user who wants to quickly get output 
without having to do any scripting.

# Installation from Github Repository

You have two options of how to do this - either is perfectly acceptable.

## Method 1:

Note for this you must have the "devtools" installed as previously recommended.

1. Open an R window and type into it: install_github("HybRIDS", username="Ward9250", ref="master")
2. Thats it!

## Method 2:

1. Click the buttons in the navigation bar called ZIP or TAR to download the files, or get them from the github repository - it is up to you whether you download the tar or zip folder and does not matter. Unzip the file.
2. Start a terminal/command-prompt and issue the command R CMD check followed by the path to where you juse unzipped the file. For example R CMD check ~/Desktop/HybRIDS. This will run a series of checks HybRIDS should pass.
3. After the check then run R CMD build --binary or R CMD build followed by the path to the repo. The first builds a binary file, the second makes a source tarball file.
4. After building the package then issuing R CMD INSTALL followed by the path to the file created from step 4.


# Installation from Binary or Source from The NBI Website

R Packages are typically distributed as either a source or binary distribution. Instructions for both are given below.

First it is essential to have a working up to date install of R on your computer. HybRIDS has been developed largely
under R version 2 but was finalised under R version 3 and used with R version 3 without issue.

You can get a stable working binary file from the ELSA website:

## On Windows

### From binary distribution

1. On Windows open the R GUI (The R icon in your programs desktop). (It should not matter if it is 32 or 64 bit).
2. Go to Packages & Data menu.
3. If you downloaded the file from the ELSA website change the dropdown menu that says "Cran Binaries" so as it selects "Local Binary Package".
4. Click the install button and then in the window that appears, select the binary file - wherever that may be on your computer. Then click the open button.

### From source distribution

1. On Windows open the R GUI (The R icon in your programs desktop). (It should not matter if it is 32 or 64 bit).
2. Go to Packages & Data menu.
3. If you downloaded the file from the ELSA website change the dropdown menu that says "Cran Binaries" so as it selects "Local source Package".
4. Click the install button and then in the window that appears, select the binary file - wherever that may be on your computer. Then click the open button.

## On a Mac

### From binary distribution

1. On mac open the R GUI (The R icon in your Mac applications).
2. Go to Packages & Data menu.
3. If you downloaded the file from the ELSA website change the dropdown menu that says "Cran Binaries" so as it selects "Local Binary Package".
4. Click the install button and then in the window that appears, select the binary file - wherever that may be on your computer. Then click the open button.

### From source distribution

1. On mac open the R GUI (The R icon in your Mac applications).
2. Go to Packages & Data menu.
3. If you downloaded the file from the ELSA website change the dropdown menu that says "Cran Binaries" so as it selects "Local source Package".
4. Click the install button and then in the window that appears, select the binary file - wherever that may be on your computer. Then click the open button.

## On Linux based systems

### From a binary distribution

1. Open R in a linux terminal and use the following command:

``` install.packages(file_name_and_path, repos = NULL, type="source") ```

Replace "file_name_and_path" with the name and location of the file.

### From binary distribution

1. Open R in a linux terminal and use the following command:

``` install.packages(file_name_and_path, repos = NULL) ```

Replace "file_name_and_path" with the name and location of the file.

# Running HybRIDS

Running HybRIDS is simple: On your machine you first need to start R. Do this by either typing R on the linux command line to start R in the terminal 
(some nicer distros like Ubuntu may have an easy to use icon/shortcut, but don't rely on that!), or if you're on Mac or Windows you will find you start 
R just like you would any other windows program or macbook app.

Once R is open, and you have installed HybRIDS as described above you simply need to invoke the library command:

``` library(HybRIDS) ```

You will be asked if you want to use the graphical user interface (GUI) for HybRIDS. If you do not need it and instead want to use HybRIDS's exported 
functions, then you can click no.
Note that saying yes to the GUI does not mean you cannot use the exported functions for use in the R console.

## Internal documentation

The documentation of each HybRIDS function can be viewed from within R by typing:
``` ?function_name_here ``` in the R console.

# Contributing

If anyone wants to contribute in the future, we welcome it, whether it is improvements over old code or new features/ideas.
Contribution can be done by the usual github cycle of sending pull requests, alternatively contact one of the authors of the paper: 

# Filing Issues

If you have trouble working the package on your data there may be many causes from a genuine bug in the code - a misreading of a file leading to downstream issues, or unsuitability of the data for the kind of analysis HybRIDS does.
HybRIDS has been developed and tested using both real datasets and simulated data. However of course that does not mean there are no datasets which may cause problems. If you run into problems, you can contact the maintainer or file an issue on the repository. 
Be descriptive and detailed in what you did so as the error can be reproduced, a sample of the data that causes the error might be needed to get to the bottom of the issue.

# Known issues:

1. Bug in block detection algorithm which finds smaller blocks at consecutively smaller thresholds at the edges of true larger blocks at higher thresholds. - We have a fix for this, but just needs to be implemented.