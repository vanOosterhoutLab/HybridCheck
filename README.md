<a name="logo"/>
<div align="center">
<a href="http://ward9250.github.io/HybridCheck" target="_blank">
<img src="http://ward9250.github.io/HybridCheck/img/Hybrid-Checklogo.png" alt="HybridCheck Logo Here"></img>
</a>
</div>

[![Gitter](https://badges.gitter.im/Join Chat.svg)](https://gitter.im/Ward9250/HybridCheck?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

[![Build Status](https://travis-ci.org/Ward9250/HybridCheck.svg?branch=master)](https://travis-ci.org/Ward9250/HybridCheck)

**An R package for simple visualisation of recombinant regions and mosaic genome structures. Including block detection and simple dating estimation.**
**Visit the website ward9250.github.io/HybridCheck**

# Introduction

The HybridCheck project is an R package that is intended to make it quick, simple and easy to script scans of recombination signal in sequence triplets. 

# Installation from Github

Note for this you must have the "devtools" installed.

1. Open an R window and type into it: `install_github("Ward9250/HybridCheck", build_vignettes=TRUE)`
2. Thats it!

# Running HybridCheck

Once R is open, and you have installed HybridCheck as described above you simply need to invoke the library command:

```R library(HybridCheck) ```

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
