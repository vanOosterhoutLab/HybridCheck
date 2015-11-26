echo " ** THIS SCRIPT WILL INSTALL THE HybridCheck R PACKAGE TO YOUR SYSTEM **"
echo " - Ensure that you have the nessecery permissions and a working R installation before you try this"
echo "   Otherwise this script will not work."

echo "Starting R..."

Rscript "Install_HybridCheck.R"

echo alias hybridcheck=\''Rscript -e "library(shiny); library(shinydashboard); runGitHub(\"Ward9250/HybridCheckApp\", launch.browser = TRUE)"'\' >> ~/.bashrc

echo " - The HybridCheck R package, and it's dependencies have been installed."
echo " - To start the HybridCheck interactive web-app, use the command 'hybridcheck' in a new terminal window."
