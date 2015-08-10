#' Reference class to store and manage triplet data.
#' @name Triplet
Triplet <-
  setRefClass("Triplet",
              
              fields = list(
                ScanAnalysis = "ANY",
                Blocks = "ANY"
              ),
              
              methods = list(
                initialize = function(){
                  
                }
              )
)