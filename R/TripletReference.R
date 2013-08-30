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
                                   },
                                 
                                 # Method for plotting the Linesplot with ggplot2 for Sequence Similarity.
                                 plotLines =
                                   function( LabelFontSize = 12, LegendFontSize = 12, Title = TRUE ) {
                                     combo <- unlist( lapply( combn( c( SequenceA, SequenceB, SequenceC ), 2, simplify=FALSE ), function(x) paste( x, collapse=":" ) ) )
                                     similarities <- as.matrix( SSTable[ , 7:9] )
                                     plotting.frame <- data.frame( basepos = rep( as.numeric( SSTable[,4] ), 3 ),
                                                                   xrange = rep( c( 1:nrow( similarities ) ) ),
                                                                   yvalues = as.vector( similarities ),
                                                                   factors = rep( 1:3, each = nrow( similarities ) ) )
                                     plot <- ggplot(plotting.frame, aes(x=basepos, y=yvalues)) + geom_line(aes(colour=factor(factors)), show_guide=T, size=0.8) +
                                       ylim(0,100) + 
                                       scale_colour_manual(name = "Pairwise Comparrisons", labels=c(combo[1], combo[2], combo[3]),values=c("yellow","purple","cyan")) +
                                       xlab("Base Position") +
                                       ylab("% Sequence Similarity")
                                     if( Title == TRUE ) {
                                       plot + 
                                         ggtitle( paste("Sequence Similarity Between Sequences for Triplet ", SequenceA, ":", SequenceB, ":", SequenceC, sep="") ) +
                                         theme( title = element_text(size = 14, colour = "black", face = "bold" ),
                                                axis.title.y = element_text( size = LabelFontSize, colour = "black" ),
                                                axis.title.x = element_text( size = LabelFontSize, colour = "black" ), 
                                                axis.text.x = element_text( size = LabelFontSize, colour = "black" ), 
                                                legend.text = element_text( size = LegendFontSize ) )
                                     } else {
                                       plot + theme( axis.title.y = element_text( size = LabelFontSize, colour = "black" ),
                                                     axis.title.x = element_text( size = LabelFontSize, colour = "black" ), 
                                                     axis.text.x = element_text( size = LabelFontSize, colour = "black" ), 
                                                     legend.text = element_text( size = LegendFontSize ) )
                                     }
                                   })
                               )