# Reference class for Triplets 
HybRIDStriplet <- setRefClass( "HybRIDStriplet",
                               fields = list( SSTableFile = "character",
                                              SSTable = function( value ) {
                                                if( missing( value ) ) {
                                                  read.table( SSTableFile )
                                                } else {
                                                  write.table( value, file = SSTableFile )
                                                }
                                              },
                                              InformativeDNALength = "numeric",
                                              SequenceA = "character",
                                              SequenceB = "character",
                                              SequenceC = "character",
                                              WindowSizeUsed = "numeric",
                                              StepSizeUsed = "numeric",
                                              SSError = "character",
                                              SSWarning = "character",
                                              Blocks = "list"
                                              ),
                               methods = list(
                                 initialize = 
                                   function( sequences ) {
                                     SSTableFile <<- tempfile( pattern = "SSTable" )
                                     SSTable <<- data.frame( WindowCenter = NA, WindowStart = NA, WindowEnd = NA,
                                                             ActualCenter = NA, ActualStart = NA, ActualEnd = NA,
                                                             AB = NA, AC = NA, BC = NA )
                                     SequenceA <<- sequences[1]
                                     SequenceB <<- sequences[2]
                                     SequenceC <<- sequences[3]
                                   })
                               )