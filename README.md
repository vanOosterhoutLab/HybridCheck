# HybirdCheck

## Introduction

HybridCheck is an R package that is intended to make it quick, simple
and easy to script scans of recombination signal in sequence triplets.

It also allows you to compute ABBA BABA tests for gene flow and some variants of
the ABBA BABA test.

## Installation instructions

Installing HybridCheck from GitHub is easy if you have "devtools" installed.
If you don't have devtools, fear not, it is easy to acquire.

```R
install.packages("devtools")
```

Once you have devtools, you can use the `install_github` function to grab and
install HybridCheck.

```R
devtools::install_github("BenJWard/HybridCheck", build_vignettes = TRUE)
```

## Using HybridCheck

Once R is open, and you have installed HybridCheck as described above you simply need to invoke the library command:

```R
library(HybridCheck)
```

### Loading sequence data

Create a new HybridCheck session by using `HC$new`. Provide a fasta file as an
argument. For example:

```R
MyAnalysis <- HC$new("~/Dropbox/MySequences.fas")
```



## Documentation

The documentation of each HybridCheck class and method function can be viewed from within R by typing:
```R ?function_name_here ``` in the R console.

New users should start with the HybridCheck User Manual vignette provided with the package in the vignettes directory.
It is also available as a PDF from the website.

To get the vignette enter in your R session:
```R
vignette("HybridCheck_user_manual")
```

A vignette on scripting with the HybridCheck library can be obtained with:
```R
vignette("Scripting_with_HybridCheck")
```

## Test Data

An example sequence alignment of 10 sequences as a DNAbin object is available in the package to use to test HybridCheck functionality.
```R
data(MySequences)
```

# HybridCheck-App

A Shiny based web-app is available to help people use this package's functionality with a GUI or on a server / cluster.
It is very new compared to this R package so feedback is always welcome and it's layout and functionality is likely to change in the future. It is located at the [following repo.](https://github.com/Ward9250/HybridCheck-App)


# Contributing

If anyone wants to contribute in the future, we welcome it, whether it is improvements over old code or new features/ideas.
Contribution can be done by the usual github cycle of sending pull requests, alternatively contact one of the authors of the paper:

# Bugs and Issues

If you have trouble working the package on your data there may be many causes from a genuine bug in the code - a misreading of a file leading to downstream issues, or unsuitability of the data for the kind of analysis HybridCheck does.
HybridCheck has been developed and tested using both real datasets and simulated data. However of course that does not mean there are no datasets which may cause problems. If you run into problems, you can contact the maintainer or file an issue on the repository.
Be descriptive and detailed in what you did so as the error can be reproduced, a sample of the data that causes the error might be needed to get to the bottom of the issue.
