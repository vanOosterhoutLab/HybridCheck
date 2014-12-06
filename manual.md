---
layout: manualsection
title: Manual Home
permalink: manual.html
manual: true
---

<embed src="http://ward9250.github.io/HybRIDS/HybRIDS_user_manual.pdf" width="500" height="375">


The HybRIDS manual
------------------

This manual is written with the intention of being your first stop for finding out how to use HybRIDS to do an
analysis of some of your sequences. Before you go any further, you should head over to our getting started page for advice 
on installing R and so HybRIDS. That information is given in a separate section and page, rather than the manual, because installing and 
getting to grips with R in detail, would be a huge and arguably unrelated segment to the manual on it's own.

The manual is divided up according to the tasks you would perform with HybRIDS, from loading the package and reading in data, to altering the options for and running the analysis, and exporting data and plotting data.

Stuck?
------

If you get stuf using R there are a great many online sources of help out there, I recommend asking questions on Stack Overflow, or searching stack overflow - your R question may have already been answered!
There is also Quick-R, R-Bloggers, the R mailing lists. R is one of those languages where googling for help (or just googling about R) is always bound to A.) help you solve your problem and B.) most likely teach you something new about R.

But what if you get stuck with a HybRIDS specific issue? Perhaps the manual is unclear in an area or you would like to do
something not covered by the manual or perhaps there is even a bug (as an open source and continually developing tool it may happen!).
In this case, and indeed even in the case where you need general R help, I encourage you to contact us by email and we'll do our best to help you.

---

An Introduction to HybRIDS's structure
--------------------------------------

If you are unfamiliar how R works as a language or not interested in how HybRIDS works you can safely ignore this, you should still be able to use it:

HybRIDS revolves around several reference classes. Reference classes are a relatively recent addition to R and violate the traditional 
R mindset as a functional language, if you're not used to programming terminology and this means nothing to you: 
The pass by reference mechanism - in a nutshell - helps us get around issues with the way R normally works as a 
functional language: In R, passing data to a function results in a copying of it once it is modified, resulting in 
large memory uses where there need not be. Note however, for many statistical tasks the origional copy mechanism is good - you don't want your origional dataset sullied and so passing a copy of the data that is changed in a function is safer.

DNA data can be large, as can the data HybRIDS generates, so passing by reference, which makes changes to the original data without copy, makes more sense and also feels more like programming in other languages like C++ or Java.

It also makes using the HybRIDS library a bit easier to use for the non R enthusiast. All the commands and steps of a HybRIDS analysis are done by calling the methods of a single HybRIDS object. So you only have one object/variable in your R workspace to worry about, and it contains everything you should need for using HybRIDS core functionality.

This manual takes you through the process of using HybRIDS in exactly such a manner - using the R command line. However for those who really do not want to touch a command line or scripting with R,
there is a Graphical user interface HybRIDS offers. This uses the gWidgets R package and the tcltk GUI toolkit in order to work. I personally recommend the scripting option, it is far more flexible and easier to use HybRIDS with other scripts and libraries. In addition, changes and new features implemented in HybRIDS will be useable immediately from the command line, whereas GUI users would have to wait a little longer until I have the chance to alter the GUI to have new buttons and etc for the new feature.

---

Finally it is worth noting that I do replicate some of the information in this manual on the Wiki Page at the HybRIDS github repository.
Information there is transferred here usually when we think it is mature and unlikely to change in the future. The wiki is then more where 
we write up information for in development stuff or debugging and other purposes, but as a result some of the info there is the same as it is 
here and in some cases it is not - always go with the instructions in this manual first as it directly applies to the master branch i.e. the stable version of the package.
