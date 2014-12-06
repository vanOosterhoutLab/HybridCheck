---
layout: manualsection
title: Testing and dating detected blocks
permalink: 08-manual.html
manual: true
published: true
status: publish
---
 

 
Having searched for blocks in the desired sequence triplets, the `BlockDating` step is executed using the `dateBlocks` method.
 
This step assigns a significance value to blocks based on the size of a block, the number of mutations observed, and the binomial probability distribution.
 
The parameters for this step are:
 
MutationRate
  : The substitution rate assumed when estimating the ages of blocks. By default it is 10e-9.
 
PValue
  : The critical alpha value when testing the significance values of blocks. By default is is 0.05.
  
BonfCorrection
  : Logical (`TRUE` / `FALSE`) when `TRUE` (default) the `PValue` will be adjusted with a Bonferroni correction.
  
DateAnyway
  : Logical (`TRUE` / `FALSE`) when `TRUE` (`FALSE` by default), all detected blocks will be kept and dated.
  : When `FALSE`, blocks that do not pass the critical alpha when testing for significance are not kept in the results
  : and are not dated.
  
  
 
 
 
 
 
