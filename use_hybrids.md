---
layout: default
title: Get and deploy HybridCheck
subtitle: Install and use HybridCheck in one of several ways.
permalink: use_HybridCheck.html
---

### You have several options for installing and using HybridCheck.

There are several different ways you can install and use HybridCheck:

1. If you don't want to use R, you may use one of the pre-packaged utilities we
   provide for Windows, OS X, and Linux.

2. We provide Docker containers on Dockerhub which provide an image, which may be
run on any machine compatible with Docker. This Docker image has R, HybridCheck and HybridCheckApp
pre-installed.

3. If you have an up to date installation of R on your machine you can install
HybridCheck straight from Github using the R devtools package.

-----

### 1. Using HybridCheck with the pre-packaged utilities for Windows 7, OSX, and Linux.

For users who want to use HybridCheck without needing to know much R, we have developed a Shiny web-app based GUI.

In addition we have provided a launcher for Windows, OS X, and Linux. So long as you have R on your machine, these launchers will make sure HybridCheck and its dependencies are installed, will fetch the current version of the Shiny files for the app from the GitHub repository, and then start up the app. You can then use HybridCheck with the GUI as detailed in the [HybridCheck User Manual.](./manual.html)

##### Operating system specific installers
<div class="container">
      <div class="row">
        <div class="col-md-4">
          <h3>Windows</h3>
          Download this installer and run to install the launcher to launch HybridCheck from the Start > All Programs menu.
          <a class="btn btn-default" href="./Installers/Windows/HybridCheckInstaller.EXE" role="button" download="install_HybridCheckLauncher.EXE">Download for Windows</a>
        </div>
        <div class="col-md-4">
          <h3>OS X</h3>
          Download this .app launcher and drag it to your Applications folder. This .app uses JavaScript for Automation and so is compatible with recent versions of OS X and has been tested on Yosemite.
          <p><a class="btn btn-default" href="./Installers/OSX/HybridCheck.app.zip" role="button" download="HybridCheckApp.zip">Download for OS X</a></p>
       </div>
        <div class="col-md-4">
          <h3>Linux</h3>
          Download this tarball and extract the folder. The folder includes two shell scripts that can be used to install HybridCheck and its dependencies, and also run HybridCheck with the webapp GUI.
          <p><a class="btn btn-default" href="./Installers/Linux/HybridCheck_Linux_Installer.zip" role="button" download="HybridCheck_Linux_Installer_Scripts.zip">Download for GNU/Linux</a></p>
        </div>
      </div>
</div>

-----

### 2. Installing and using a Docker container.

[Docker](https://www.docker.com) is an open platform for building, shipping and
running distributed applications. We provide docker containers with R, HybridCheck,
HybridCheckApp and any dependencies installed. These containers can be installed from
[DockerHub](https://hub.docker.com) using the docker utility.

We provide two containers:

#### 1. [ward9250/hybridcheck](https://hub.docker.com/r/ward9250/hybridcheck/)
A Docker container with R, HybridCheck and it's dependencies installed. Use this
if you just want to use the HybridCheck package from the R console.


#### 2. [ward9250/dockerized-hybridcheckapp](https://hub.docker.com/r/ward9250/dockerized-hybridcheckapp/)
A Docker container with R and HybridCheck installed, and which also has Shiny
and HybridCheckApp installed. Use this container if you want to use the web-app
designed to provide a Graphical User Interface for the R package.

-----

### 3. Installing from the R console.

If an installer for your operating system is not available or working for you,
which may happen due to permissions. You can install HybridCheck and every one
of it's dependencies to you R library manually by pasting the following commands
into an R console:

```R
chooseCRANmirror(ind = 83)
pkg <- c('devtools', 'shiny', 'ggplot2', 'grid', 'gridExtra', 'ape', 'png')
new.pkg <- pkg[!(pkg %in% installed.packages())]
if(length(new.pkg)){install.packages(new.pkg)}
source('http://bioconductor.org/biocLite.R')
if(!('Biostrings' %in% installed.packages())){ biocLite('Biostrings', ask = F)}
if(!'IRanges' %in% installed.packages()){biocLite('IRanges', ask = F)}
if(!'HybridCheck' %in% installed.packages()){library(devtools); install_github('Ward9250/HybridCheck', ref = 'master')}
```
This bit of code explicitly downloads the dependencies HybridCheck needs in
addition to the HybridCheck package itself.

You can then run HybridCheck with the web-app by pasting the following in the R console:

```R
library(shiny)
runGitHub('HybridCheckApp', 'Ward9250', launch.browser = TRUE)
```

-----

