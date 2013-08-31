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
                                              FullDNALength = "numeric",
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
                                   function( sequences, fullseqlength ) {
                                     SSTableFile <<- tempfile( pattern = "SSTable" )
                                     SSTable <<- data.frame( WindowCenter = NA, WindowStart = NA, WindowEnd = NA,
                                                             ActualCenter = NA, ActualStart = NA, ActualEnd = NA,
                                                             AB = NA, AC = NA, BC = NA )
                                     SequenceA <<- sequences[1]
                                     SequenceB <<- sequences[2]
                                     SequenceC <<- sequences[3]
                                     FullDNALength <<- fullseqlength
                                   },
                                 
                                 # Method for plotting the Linesplot with ggplot2 for Sequence Similarity.
                                 plotLines =
                                   function( LabelFontSize = 12, LegendFontSize = 12, title = TRUE ) {
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
                                     if( title == TRUE ) {
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
                                   },
                                 
                                 #Plotting method for the rainbow bars in ggplot2
                                 plotBars =
                                   function( mosaicscale = 500, pyramidleg = TRUE, labfontsize = 12, legfontsize = 12, xlabel = TRUE, title = TRUE ) {
                                     # Now let's generate the reference colour palettes.
                                     RefA <- expand.grid(contigb = seq(0, 100, by = 1), contigc = seq(0, 100, by = 1))
                                     RefA <- within(RefA, mix <- rgb(green = contigb, red = 100, blue = contigc, maxColorValue = 100))
                                     RefB <- expand.grid(contiga = seq(0, 100, by = 1), contigc = seq(0, 100, by = 1))
                                     RefB <- within(RefB, mix <- rgb(green = 100, red = contiga, blue = contigc, maxColorValue = 100))
                                     RefC <- expand.grid(contiga = seq(0, 100, by = 1), contigb = seq(0, 100, by = 1))
                                     RefC <- within(RefC, mix <- rgb(green = contigb, red = contiga, blue = 100, maxColorValue = 100))
                                     # Now figure out the scale and data to go into each vertical bar:
                                     div <- FullDNALength / mosaicscale
                                     bpstart <- seq(from = 1, to = FullDNALength, by = div)
                                     bpend <- seq(from=div, to = FullDNALength, by = div)
                                     bpX <- round( bpstart +  ( div / 2 ) )
                                     frame <- data.frame(bpstart = bpstart, bpend = bpend, bpX = bpX )
                                     rm(bpstart, bpend)
                                     # Now we go through each vertical bar, and for each one, we find the SSvalues to go in there, and we average them.
                                     ss <- SSTable
                                     AB <- round( apply( frame, 1, function(x) vertbar_create( ss, x, 7 ) ) )
                                     AC <- round( apply( frame, 1, function(x) vertbar_create( ss, x, 8 ) ) )
                                     BC <- round( apply( frame, 1, function(x) vertbar_create( ss, x, 9 ) ) )
                                     rm( ss )
                                     frame$AB <- AB
                                     frame$AC <- AC
                                     frame$BC <- BC
                                     frame$X <- 1:nrow(frame)
                                     rm( AB, AC, BC )
                                     frame$A_mix <- apply( frame, 1, function(x) col_deter( c( as.numeric(x[4]), as.numeric(x[5]) ), RefA ) )
                                     frame$B_mix <- apply( frame, 1, function(x) col_deter( c( as.numeric(x[4]), as.numeric(x[6]) ), RefB ) )
                                     frame$C_mix <- apply( frame, 1, function(x) col_deter( c( as.numeric(x[5]), as.numeric(x[6]) ), RefC ) )
                                     plottingFrame <- data.frame( X = frame$X, Y = rep( c(3, 2, 1), each = mosaicscale), colour = c(frame$A_mix, frame$B_mix, frame$C_mix))
                                     bars <- ggplot(plottingFrame, aes(x = X, y = as.factor(Y)) ) +
                                       geom_raster( aes( fill = colour ) ) + scale_fill_identity() +
                                       xlab("Approximate Base Position") +
                                       ylab( "Sequence Name" ) +
                                       scale_x_continuous( breaks = c(seq( from = 1, to = mosaicscale, by = mosaicscale / 10 ), mosaicscale), labels = c(bpX[seq( from = 1, to = mosaicscale, by = mosaicscale / 10 )], max(bpX)) ) + 
                                       scale_y_discrete( labels = c(Sequence3, Sequence2, Sequence1))
                                     if( xlabel == TRUE ){
                                       bars <- bars + theme( title = element_text(size = 14, colour = "black", face = "bold" ),
                                                     axis.text.y = element_text( colour="black", size = 12 ),
                                                     axis.title.y = element_text( size = labfontsize, colour = "black" ),
                                                     axis.text.x = element_text( colour="black", size = 12),
                                                     axis.title.x = element_text( colour = "black", size = labfontsize)
                                       )
                                     } else {
                                       bars <- bars + theme( title = element_text(size = 14, colour = "black", face = "bold" ),
                                                     axis.text.y = element_text( colour="black", size = 12 ),
                                                     axis.title.y = element_text( size = labfontsize, colour = "black" ),
                                                     axis.text.x = element_blank(),
                                                     axis.title.x = element_blank()
                                       )
                                     }
                                     if( title == T ) {
                                       bars <- bars + ggtitle( paste("Sequence Similarity Between Sequences for Triplet ", SequenceA, ":", SequenceB, ":", SequenceC, sep="") )
                                     }
                                     if( pyramidleg == T ) {
                                       legend <- readPNG( system.file( "extdata/rgblegend.png", package="HybRIDS" ), TRUE )
                                       if ( names( dev.cur() ) == "windows" ) {
                                         # windows device doesn’t support semi-transparency so we’ll need
                                         # to flatten the image
                                         transparent <- legend[,,4] == 0
                                         legend <- as.raster( legend[,,1:3] )
                                         legend[transparent] <- NA
                                       }
                                       legendgrob <- rasterGrob( image=legend )
                                       bars <- arrangeGrob( bars, legendgrob, widths = c( 1, 0.13 ), ncol=2)
                                       print( bars )
                                     } else {
                                       print( bars )
                                     }
                                   },
                                 
                                 # Method for putative block detection.
                                 putativeBlockFind = 
                                   function( autodetect = TRUE, sd.stringency = 2, manual.thresholds = c( 90 ), manual.fallback = T ) {
                                     if( autodetect == TRUE ) {
                                       cat( "Using the autodetect thresholds method...\n" )
                                       cat( "Deciding on suitable thresholds...\n" )
                                       Thresholds <- autodetect.thresholds( SSTable, sd.stringency, manual.thresholds, manual.fallback )
                                       # Results in a list of thresholds for AB, AC and BC.
                                     } else {
                                       Thresholds <- list( manual.thresholds, manual.thresholds, manual.thresholds )
                                     }
                                     names(Thresholds) <- c( paste( SequenceA, SequenceB ,sep=":" ), paste( SequenceA, SequenceC, sep=":" ), paste( SequenceB, SequenceC, sep=":" ) )
                                     cat( "Now beginning Block Search...\n\n" )
                                     Blocks <<- lapply( 1:3, function(i) block.find( SSTable[,c( 1:6, 6+i )], Thresholds[[i]] ) )
                                     names(Blocks) <<- names(Thresholds) <- c( paste( SequenceA , SequenceB, sep=":" ), paste( SequenceA, SequenceC, sep=":" ), paste( SequenceB, SequenceC, sep=":" ) )
                                     
                                   }
                                 )
                               )