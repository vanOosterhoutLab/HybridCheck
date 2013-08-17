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

Installing the bleeding edge versions from Github is slightly different to installing the source or binary tarballs resent on the NBI website - which are not updated as often.
The Github repository is updated and watched almost daily by the author and so is recommended. The process is the same on all machines and operating systems.
The main difference is that the Github repo is raw source code, unpackaged into a source tarball or a binary file. It is cloned to a PC as directory/set of files that follow the R layout for packages.
Therefore there is an added step that once you've "cloned this repo" to your computer i.e. downloaded the files, there are a few steps extra just to put it into an R source tarball or build it before you can then install it in the usual manner.

1. First make sure that R is installed correctly for your system, here I must direct you to the R website and the CRAN. If youknow how to install 'normal' programs on your PC or Mac this should be trivial.
2. Next click the button on the site to clone the repo to your computer or if you have git installed and are comfortable using it, you can issue a command-line call to clone to your computer using git.
3. Start a terminal/command-prompt and issue the command ``` R CMD check ``` followed by the path to where you downloaded the repo. For example ``` R CMD check ~/Desktop/HybRIDS ```. This will run a series of checks HybRIDS should pass.
4. After the check then run ``` R CMD build --binary ``` or ``` R CMD build ``` followed by the path to the repo. The first builds a binary file, the second makes a source tarball file.
5. After building the package then issuing ``` R CMD INSTALL ``` followed by the path to the file created from step 4.

Currently it does not really matter if you install HybRIDS on your R system from binary or source since there is - at present - no object files or C/FORTRAN code in the package, it is pure R.
This may change in the future however.


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

## The exported functions

The documentation of each can be viewed from within R by typing:
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