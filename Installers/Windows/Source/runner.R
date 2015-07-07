message("Checking all required packages installed.")

chooseCRANmirror(ind = 83)

pkg <- c('devtools', 'shiny', 'shinydashboard', 'ggplot2', 'grid', 'gridExtra', 'ape', 'png')
new.pkg <- pkg[!(pkg %in% installed.packages())]
if(length(new.pkg)){install.packages(new.pkg)}
source('http://bioconductor.org/biocLite.R')
if(!('Biostrings' %in% installed.packages())){biocLite('Biostrings', ask = F)}
if(!'IRanges' %in% installed.packages()){biocLite('IRanges', ask = F)}
if(!'HybridCheck' %in% installed.packages()){library(devtools); install_github('Ward9250/HybridCheck', ref = 'master')}
library(shiny)
runGitHub('HybridCheckApp', 'Ward9250', launch.browser = TRUE)
