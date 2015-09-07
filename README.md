<a name="logo"/>
<div align="center">
<a href="http://ward9250.github.io/HybridCheck">
<img src="http://ward9250.github.io/HybridCheck/img/HybridCheckLogo.png" height="250" alt="HybridCheck Logo Here"></img>
</a>
</div>
<p align="center">
    <a href="https://gitter.im/Ward9250/HybridCheck?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge">
        <img src="https://badges.gitter.im/Join Chat.svg" alt="Gitter">
    </a>
    <a href="https://travis-ci.org/Ward9250/HybridCheck">
        <img src="https://travis-ci.org/Ward9250/HybridCheck.svg?branch=master" alt="Build Status">
    </a>
    <a href="https://waffle.io/Ward9250/HybridCheck">
        <img src="https://badge.waffle.io/Ward9250/HybridCheck.png?label=ready&title=Ready" alt="Stories in Ready">
    </a>
    <a href="https://waffle.io/Ward9250/HybridCheck">
        <img src="https://badge.waffle.io/Ward9250/HybridCheck.png?label=In%20Progress&title=In%20Progress" alt="Stories in Ready">
    </a>
</p>

# File format converters

The HybridCheck project is an R package that is intended to make it quick,
simple and easy to script scans of recombination signal in sequence triplets.

When HybridCheck was written, it was designed in response to specific needs of
several of my (Ward9250) academic projects. These projects involved analyzing
large assembled contigs or large consensus sequences, generated from reads
aligned to references.

As more people became interested and asked about HybridCheck here at the NRP,
we learned that not everyone has FASTA files like that available for all their organisms/isolates/individuals/etc.

Sometimes it was easy to convert or process their data into an aligned FASTA
format. Sometimes it was more difficult.

I (Ward9250) am working on a second version of HybridCheck. This will add more
features, fix stuff in the codebase I hated about the first version, refactor
code, and improve performance. There is no set ETA for version 2 yet, because of
other projects and PhD commitments like thesis writing.

This second version shall be able to parse and make use of multiple file formats
and convert between some formats. Especially those formats dedicated to variant
and mutation data.

This repository forms an intermediate solution for converting file formats to a
form suitable for HybridCheck. Scripts will be added over time as they are created.

# Scripts

| Script Name    | Description                                                  |
|----------------|--------------------------------------------------------------|
| MAF_2_FASTA.py | A Python script, which converts Multiple Alignment format    |
|                | (MAF), into FASTA files suitable for parsing by HybridCheck. |
|                | At the time of writing, this script depends on an up to      |
|                | date BioPython installation with functionality to read in    |
|                | MAF files. This functionality is provided by the             |
|                | Bio.AlignIO.MafIO module.                                    |


# Bugs and Issues

HybridCheck has and associated files, have been developed and tested using both
real datasets and simulated data.
However of course that does not mean there are no datasets which may cause problems.
If you run into problems, you can contact the maintainer or file an issue on the repository.
Be descriptive and detailed in what you did so as the error can be reproduced, a sample of the data that causes the error might be needed to get to the bottom of the issue.

# Contribution

We welcome contribution to the project in all forms, contact the maintainers
listed on the website, or drop a line on the HybridCheck Gitter chatroom.
