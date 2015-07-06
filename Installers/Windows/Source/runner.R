message("Checking all required packages installed.")

chooseCRANmirror(ind = 83)

if(!"devtools" %in% rownames(installed.packages())){
  install.packages("devtools")
}
if(!"HybridCheck" %in% rownames(installed.packages())){
  library(devtools)
  install_github("Ward9250/HybridCheck")
}
if(!"shiny" %in% rownames(installed.packages())){
  install.packages("shiny")
}
if(!"shinyBS" %in% rownames(installed.packages())){
  library(devtools)
  install_github("ebailey78/shinyBS")
}

library(shiny)
library(shinyBS)
runGitHub("HybridCheckApp", "Ward9250", launch.browser = TRUE)
