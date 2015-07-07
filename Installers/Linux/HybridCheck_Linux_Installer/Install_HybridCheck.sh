Rscript "Install_HybridCheck.R"
echo alias "hybridcheck=\"Rscript -e \"library(shiny); library(shinydashboard); runGitHub('HybridCheckApp', 'Ward9250', launch.browser = TRUE)\"\"" >> ~/.bashrc
echo " - The HybridCheck R package, and it's dependencies have been installed."
echo " - To start the HybridCheck interactive web-app, use the command 'hybridcheck' on the command line."
