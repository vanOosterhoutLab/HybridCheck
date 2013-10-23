.onAttach <- function(...) {
  cat("\nTACTTTGTACCTAAGTATGCATTACGTTACGTTAGTAGCTGGACCTAGTAAATCGGA     
,--.  ,--.         ,--.   ,------. ,--.,------.   ,---.
|  '--'  |,--. ,--.|  |-. |  .--. '|  ||  .-.  \\ '   .-'
|  .--.  | \\  '  / | .-. '|  '--'.'|  ||  |  \\  :`.  `-.
|  |  |  |  \\   '  | `-' ||  |\\  \\ |  ||  '--'  /.-'    |
`--'  `--'.-'  /    `---' `--' '--'`--'`-------' `-----'
          `---'
ATGAAACATGGATTCATACG - Version 1.0 - CGACCTGGATCATTTAGCCT\n\n
Hybridisation, Recombination and Introgression Detection
                   and Dating Package.\n
-----------------=========****=========------------------
Cite: TBD
Licence: GPL (Like R and most packages).
http://ward9250.github.io/HybRIDS/
-----------------=========****=========------------------\n")
  cat("\nEnter HybRIDS_GUI() at the console prompt for the interactive GUI.\n\n")
}

#' Start a GUI for HybRIDS
#' 
#' @export
HybRIDS_GUI <- function(toolkit="tcltk"){
  require(gWidgets)
  require(gWidgetstcltk)
  options(guiToolkit=toolkit)
  invisible(HybRIDS_gui_generator$new())
}