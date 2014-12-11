---
layout: default
title: Get and deploy HybRIDS
subtitle: Install and use HybRIDS in one of several ways.
permalink: use_hybrids.html
---
## Use HybRIDS in R scripts
HybRIDS is an open source package of R code, and so if you are familiar with R you can install the package in several ways and you can get programming ans scripting R code with it right away. For example you could install the latest version with devtools:

```R
install.packages("devtools")
library(devtools)
install_github("Ward9250/HybRIDS", build_vignettes=TRUE)
library(HybRIDS)
```

-----

## Using HybRIDS with a Graphical User Interface

### Deploying the Shiny app on a (web)server
However, for people who want to use it without needing to know much R, we have developed a Shiny web-app based GUI, this could be deployed on a server and logged into by many people in a lab for instance, or it can be used locally on a machine using ShinyApps. For more information about deploying shiny apps, visit the [website.](http://shiny.rstudio.com)

-----

### Using the App on your local machine with the help of a launcher.
We have provided a launcher for Windows and OSX. So long as you have R on your machine, these launchers will make sure HybRIDS and it's dependencies are installed, will fetch the current version of the Shiny files for the GUI from the GitHub repository, and then start up the app. You can then use HybRIDS with the GUI as detailed in the [HybRIDS User Manual.]()

<div class="container">
      <!-- Example row of columns -->
      <div class="row">
        <div class="col-md-4">
          <h3>Windows</h3>
          Download this installer and run to install the launcher to launch HybRIDS from the Start > All Programs menu.
          <a class="btn btn-default" href="./Windows/LaunchHybRIDS_installer.EXE" role="button" download="install_HybRIDSLauncher.EXE">Download for Windows</a>
        </div>
        <div class="col-md-4">
          <h3>OSX</h3>
          Download this .app launcher and drag it to your Applications folder.
          <p><a class="btn btn-default" href="#" role="button">Download for OSX</a></p>
       </div>
        <div class="col-md-4">
          <h3>Linux</h3>
          Donec sed odio dui. Cras justo odio, dapibus ac facilisis in, egestas eget quam. Vestibulum id ligula porta felis euismod semper. Fusce dapibus, tellus ac cursus commodo, tortor mauris condimentum nibh, ut fermentum massa justo sit amet risus.
          <p><a class="btn btn-default" href="#" role="button">View details &raquo;</a></p>
        </div>
      </div>
</div>

