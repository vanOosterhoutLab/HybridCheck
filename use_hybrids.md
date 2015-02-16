---
layout: default
title: Get and deploy HybRIDS
subtitle: Install and use HybRIDS in one of several ways.
permalink: use_hybrids.html
---

## Using HybRIDS with a Graphical User Interface

For users who want to use HybRIDS without needing to know much R, we have developed a Shiny web-app based GUI.

In addition we have provided a launcher for Windows, OS X, and Linux. So long as you have R on your machine, these launchers will make sure HybRIDS and its dependencies are installed, will fetch the current version of the Shiny files for the app from the GitHub repository, and then start up the app. You can then use HybRIDS with the GUI as detailed in the [HybRIDS User Manual.](./manual.html)

### Deploying the Shiny app on a (web)server
In addition to being run on a local machine, the HybRIDS app could be deployed on a server and logged into by many people in a lab for instance. For more information about deploying shiny apps, visit the [website.](http://shiny.rstudio.com)

## Use HybRIDS in R scripts
HybRIDS is an open source package of R code, and so if you are familiar with R you can install the package in several ways and you can get programming ans scripting R code with it right away. For example you could install the latest version with devtools:

```R
install.packages("devtools")
library(devtools)
install_github("Ward9250/HybRIDS", build_vignettes=TRUE)
library(HybRIDS)
```

-----

<div class="container">
      <div class="row">
        <div class="col-md-4">
          <h3>Windows</h3>
          Download this installer and run to install the launcher to launch HybRIDS from the Start > All Programs menu.
          <a class="btn btn-default" href="./Windows/LaunchHybRIDS_installer.EXE" role="button" download="install_HybRIDSLauncher.EXE">Download for Windows</a>
        </div>
        <div class="col-md-4">
          <h3>OS X</h3>
          Download this .app launcher and drag it to your Applications folder. This .app uses JavaScript for Automation and so is compatible with recent versions of OS X and has been tested on Yosemite. 
          <p><a class="btn btn-default" href="./OSX/HybRIDSapp.zip" role="button" download="HybRIDS.app">Download for OS X</a></p>
       </div>
        <div class="col-md-4">
          <h3>Linux</h3>
          Download this tarball and extract the folder. The folder includes two shell scripts that can be used to install HybRIDS and its dependencies, and also run HybRIDS with the webapp GUI. 
          <p><a class="btn btn-default" href="./Linux/HybRIDS_Launcher.tar.gz" role="button" download="HybRIDS_Launcher.tar.gz">Download for GNU/Linux</a></p>
        </div>
      </div>
</div>

