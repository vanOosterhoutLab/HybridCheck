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

## Using HybRIDS with a Graphical User Interface
However For people who want to use it without needing to know much R, we have developed a Shiny web-app based GUI, this could be deployed  

<div class="container">
      <!-- Example row of columns -->
      <div class="row">
        <div class="col-md-4">
          <h3>Windows</h3>
          Donec id elit non mi porta gravida at eget metus. Fusce dapibus, tellus ac cursus commodo, tortor mauris condimentum nibh, ut fermentum massa justo sit amet risus. Etiam porta sem malesuada magna mollis euismod. Donec sed odio dui.
          <p><a class="btn btn-default" href="#" role="button">View details &raquo;</a></p>
        </div>
        <div class="col-md-4">
          <h3>OSX</h3>
          Donec id elit non mi porta gravida at eget metus. Fusce dapibus, tellus ac cursus commodo, tortor mauris condimentum nibh, ut fermentum massa justo sit amet risus. Etiam porta sem malesuada magna mollis euismod. Donec sed odio dui.
          <p><a class="btn btn-default" href="#" role="button">View details &raquo;</a></p>
       </div>
        <div class="col-md-4">
          <h3>Linux</h3>
          Donec sed odio dui. Cras justo odio, dapibus ac facilisis in, egestas eget quam. Vestibulum id ligula porta felis euismod semper. Fusce dapibus, tellus ac cursus commodo, tortor mauris condimentum nibh, ut fermentum massa justo sit amet risus.
          <p><a class="btn btn-default" href="#" role="button">View details &raquo;</a></p>
        </div>
      </div>
</div>

