SSAnalysisSettings <- setRefClass("SSAnalysisSettings",
                                  
                                  fields = list(
                                    WindowSize = "integer",
                                    StepSize = "integer"
                                    ),
                                  
                                  methods = list(
                                    initialize =
                                      function(){
                                        WindowSize <<- 100L
                                        StepSize <<- 1L
                                      },
                                    
                                    getWindowSize =
                                      function(){
                                        return(WindowSize)
                                      },
                                    
                                    getStepSize =
                                      function(){
                                        return(StepSize)
                                      },
                                    
                                    setWindowSize =
                                      function(value){
                                        if(!is.integer(value)){stop("Error: Provide one integer value as a Window Size.")}
                                        WindowSize <<- value
                                      },
                                    
                                    setStepSize =
                                      function(value){
                                        if(length(value) != 1 || !is.integer(value)){stop("Error: Provide one integer value as a Step Size.")}
                                        StepSize <<- value
                                      },
                                    
                                    setSettings =
                                      function(...){
                                        settings <- list(...)
                                        parameters <- names(settings)
                                        for(i in 1:length(settings)){
                                          if(parameters[i] == "WindowSize"){
                                            setWindowSize(settings[[i]])
                                          }
                                          if(parameters[i] == "StepSize"){
                                            setStepSize(settings[[i]])
                                          }
                                        }
                                      },
                                    
                                    showSettings =
                                      function(){
                                        message("Settings Sequence Scans of Triplets")
                                        return(paste('Size of the sliding window in base pairs (WindowSize): ', getWindowSize(),
                                                     '\n\nNumber of base pairs to move the sliding window on\n\teach iteration of the scan (StepSize): ', getStepSize(), sep=""))
                                      }
                                    )
)