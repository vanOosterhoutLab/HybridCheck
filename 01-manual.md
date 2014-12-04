---
layout: default
title: Getting Started
subtitle: Get R and HybRIDS up and running on your machine.
permalink: 01-manual.html
manual: true
---

Introduction
------------
					
HybRIDS is written as an R package and so in order to use it you need a few things:
First you need the core R Environment installed on your machine. The reason being if you want to write high level language like R, Python, Perl or Ruby,
The software required to interpret the language needs to be installed.
You then need the required R packages that the HybRIDS package relies on or gains functionality when used in tandem with HybRIDS.
					
***

Installing R
------------
					
R is simple to install for most operating systems. If you use a Mac or PC or Linux system, the process for installing it is like for installing any other software on that system. 
					
					
First head over to [The Comprehensive R Archive Network [CRAN]](http://cran.r-project.org/index.html). The CRAN maintains R and officially supported packages.
This is where you will get the version of R that is appropriate for your system. At the top of the page choose one of the links "Download R for XYZ" where XYZ is Windows, Linux or Mac OSX - depending on your system.
The next page will display downloads for different versions of R - at the time of writing you will want either version 3.X.X or 2.15.X.
					
					
The 3.X versions of R mark a major release landmark and not all packages have been updated for 3.X.X just yet.
The packages HybRIDS depends on - and HybRIDS itself, can run on either the 3.X.X series or 2.15.X.
If you are on a Windows system pick a version and it will download you a .exe file you just need to run and follow the instructions and it will install on your system. On the Mac it will do similarly.
Simply go through the installer clicking 'next' and keep the defaults if you are unsure on what to do. But really it is just like installing any other program on a Windows of Mac.
If you use a Linux system is is a little more complex but still familiar to the Linux user who has installed software on it before. Instructions for the main distros are on the site and in a nutshell you must use the package manager for your distro to download and keep R up to date. If this is not possible then it will be necessary to compile R from source - in the case it is best to console the CRAN site and R community forums should you get stuck.

***
								
Installing HybRIDS
------------------
					
If you are on a Windows machine it is highly recommended that you install RTools. It can be found [here](http://cran.r-project.org/index.html).
					
There are a few different ways you might install HybRIDS to your R library. We recommend installation from GitHub is recommended because it is the easiest and can be done with only 3 commands to the R console.
				

## Installing from the GitHub Repository

This method should work no matter what OS you are using.

To install HybRIDS from the GitHub repository executing the following code - snippet into an R session should get HybRIDS installed with minimum difficulty:

```R
install.packages("devtools", dep=T)
library(devtools)
install_github("Ward9250/HybRIDS")
```

This installs a package by R community contributor Hadley Wickham which allows you to install R packages from their GitHub Repositories with minimum fuss. Hadley Wickham has contributed much to the R community and I thoroughly recomend you check out his [GitHub repository](https://github.com/hadley).
					
Now you have the latest version of HybRIDS installed in your R library. By default the above command install the master branch of the HybRIDS repository. This is the "safest" option by which we mean it is the version we are by far most confident in to contain no errors or bugs. Anything new that is written for HybRIDS first is written in it's own branch, before being pulled to the <b>devel</b> branch. When we are then sure this is stable it is merged into the master branch.

You may wish to use some of the devel features before they are merged into the master branch. In this case you would replace the last line of code in the box above with the following:

```R
install_github("Ward9250/HybRIDS", ref="devel")
```

***

Starting/Using HybRIDS
----------------------

Once you have installed HybRIDS as described above you simply need to invoke the library command as for any other installed R package:

R```
library(HybRIDS)
```