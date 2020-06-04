# HybridCheck

## Introduction

HybridCheck is an R package that is intended to make it quick, simple
and easy to script analyses scanning for recombination signal between haplotype
sequences.

It can plot such scans:

![hybridcheck-plot](exampleplot.jpeg)

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
devtools::install_github("vanOosterhoutLab/HybridCheck", build_vignettes = TRUE)
```

## Documentation

### User manual

HybridCheck comes with a user manual vignette, that you can access from an R session.

Enter:

```R
vignette("HybridCheck-manual")
```

and hit return. A new window should pop up with the manual inside.
It will guide you through loading your sequence data and running an analysis.

### API reference documentation

The documentation of each HybridCheck class and method/function can be viewed from within R by typing:

```R
?function_name_here
```

in the R console.

## Test dataset

An example sequence alignment of 10 sequences as a DNAbin object is available in the package to use to test HybridCheck functionality. For more details, see the user manual.

```R
data(MySequences)
```

## Contact

If you would like to contribute to HybridCheck, or want to ask a question or report a bug.
You can either use GitHub Issues and Pull Requests in the usual way, or you can contact either
of the two maintainers below:

- Dr. Ben J. Ward <benjward@protonmail.com> 
- Prof. Cock van Oosterhout <c.van-oosterhout@uea.ac.uk>

