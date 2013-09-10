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
                                              ContigNames = "character",
                                              WindowSizeUsed = "numeric",
                                              StepSizeUsed = "numeric",
                                              SSError = "character",
                                              SSWarning = "character",
                                              BlocksWarning = "character",
                                              Blocks = "list"
                                              ),
                               methods = list(
                                 initialize = 
                                   function( sequences, fullseqlength ) {
                                     SSTableFile <<- tempfile( pattern = "SSTable" )
                                     SSTable <<- data.frame( WindowCenter = NA, WindowStart = NA, WindowEnd = NA,
                                                             ActualCenter = NA, ActualStart = NA, ActualEnd = NA,
                                                             AB = NA, AC = NA, BC = NA )
                                     ContigNames <<- c( sequences[1], sequences[2], sequences[3] )
                                     FullDNALength <<- fullseqlength
                                   },
                                 
                                 # Method for plotting the Linesplot with ggplot2 for Sequence Similarity.
                                 plotLines =
                                   function( ... ) {
                                     Parameters <- list( ... )
                                     if( "Title" %in% Parameters ) {
                                       if(!is.logical(Parameters$Title)) stop("The Title parameter must be logical/boolean.")
                                       Title <- Parameters$Title
                                     } else {
                                       Title <- TRUE
                                     }
                                     if( "LabelFontSize" %in% names(Parameters)) {
                                       if(!is.numeric(Parameters$LabelFontSize)) stop("Parameter LabelFontSize must be numeric")
                                       LabelFontSize <- Parameters$LabelFontSize
                                     } else {
                                       LabelFontSize <- 12
                                     }
                                     if( "LegendFontSize" %in% names(Parameters)) {
                                       if(!is.numeric(Parameters$LegendFontSize)) stop("Parameter LegendFontSize must be numeric")
                                       LegendFontSize <- Parameters$LegendFontSize
                                     } else {
                                       LegendFontSize <- 12
                                     }
                                     combo <- unlist( lapply( combn( ContigNames, 2, simplify=FALSE ), function(x) paste( x, collapse=":" ) ) )
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
                                       plot <- plot + 
                                         ggtitle( paste("Sequence Similarity Between Sequences for Triplet ", ContigNames[1], ":", ContigNames[2], ":", ContigNames[3], sep="") ) +
                                         theme( title = element_text(size = 14, colour = "black", face = "bold" ),
                                                axis.title.y = element_text( size = LabelFontSize, colour = "black" ),
                                                axis.title.x = element_text( size = LabelFontSize, colour = "black" ), 
                                                axis.text.x = element_text( size = LabelFontSize, colour = "black" ), 
                                                legend.text = element_text( size = LegendFontSize ) )
                                     } else {
                                       plot <- plot + theme( axis.title.y = element_text( size = LabelFontSize, colour = "black" ),
                                                     axis.title.x = element_text( size = LabelFontSize, colour = "black" ), 
                                                     axis.text.x = element_text( size = LabelFontSize, colour = "black" ), 
                                                     legend.text = element_text( size = LegendFontSize ) )
                                     }
                                     return( plot )
                                   },
                                 
                                 #Plotting method for the rainbow bars in ggplot2
                                 plotBars =
                                   function( exportDat = FALSE, ... ) {
                                     Parameters <- list( ... )
                                     if( "Mosaic.Scale" %in% names(Parameters) ) {
                                       if(!is.integer(Parameters$Mosaic.Scale) || length(Parameters$Mosaic.Scale)) stop("The Mosaic.Scale Parameter must be an integer (e.g. 500L)")
                                       Mosaic.Scale <- Parameters$Mosaic.Scale
                                     } else {
                                       Mosaic.Scale <- 500L
                                     }
                                     if( "Legend" %in% names(Parameters) ) {
                                       if( !is.logical( Parameters$Legend) || length(Parameters$Legend) > 1) stop("The Legend plot parameter must be a boolean value.")
                                       Legend <- Parameters$Legend
                                     } else {
                                       Legend <- TRUE
                                     }
                                     if( "LabelFontSize" %in% names(Parameters) ) {
                                       if( !is.numeric(Parameters$LabelFontSize) || length(Parameters$LabelFontSize) > 1) stop("LabelFontSize parameter must be a numerical value.")
                                       LabelFontSize <- Parameters$LabelFontSize
                                     } else {
                                       LabelFontSize <- 12
                                     }
                                     if( "LabelBP" %in% names(Parameters) ) {
                                       if( !is.logical(Parameters$LabelBP) || length(Parameters$LabelBP) > 1) stop("LabelBP parameter must be a boolean value.")
                                       LabelBP <- Parameters$LabelBP
                                     } else {
                                       LabelBP <- TRUE
                                     }
                                     if( "Title" %in% names(Parameters) ) {
                                       if( !is.logical(Parameters$Title) || length(Parameters$Title) > 1) stop("Title parameter must be a boolean value.")
                                       Title <- Parameters$Title
                                     } else {
                                       Title <- TRUE
                                     }
                                     if( "TickSize" %in% names(Parameters) ) {
                                       if( !is.numeric(Parameters$TickSize) || length(Parameters$TickSize) > 1 ) stop("Title TickSize must be a numerical value.")
                                       TickSize <- Parameters$TickSize
                                     } else {
                                       TickSize <- 12
                                     }
                                     if( "TickColour" %in% names(Parameters) ) {
                                       if( !is.character(Parameters$TickColour) || length(Parameters$TickColour) > 1 ) stop("Title TickColour must be a numerical value.")
                                       TickColour <- Parameters$TickColour
                                     } else {
                                       TickColour <- "black"
                                     }
                                     # Now let's generate the reference colour palettes.
                                     RefA <- expand.grid(contigb = seq(0, 100, by = 1), contigc = seq(0, 100, by = 1))
                                     RefA <- within(RefA, mix <- rgb(green = contigb, red = 100, blue = contigc, maxColorValue = 100))
                                     RefB <- expand.grid(contiga = seq(0, 100, by = 1), contigc = seq(0, 100, by = 1))
                                     RefB <- within(RefB, mix <- rgb(green = 100, red = contiga, blue = contigc, maxColorValue = 100))
                                     RefC <- expand.grid(contiga = seq(0, 100, by = 1), contigb = seq(0, 100, by = 1))
                                     RefC <- within(RefC, mix <- rgb(green = contigb, red = contiga, blue = 100, maxColorValue = 100))
                                     # Now figure out the scale and data to go into each vertical bar:
                                     div <- FullDNALength / Mosaic.Scale
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
                                     plottingFrame <- data.frame( X = frame$X, Y = rep( c(3, 2, 1), each = Mosaic.Scale), colour = c(frame$A_mix, frame$B_mix, frame$C_mix))
                                     if( exportDat == T ) {
                                       return(frame)
                                     } else {
                                       bars <- ggplot(plottingFrame, aes(x = X, y = as.factor(Y)) ) +
                                         geom_raster( aes( fill = colour ) ) + scale_fill_identity() +
                                         xlab("Approximate Base Position") +
                                         ylab( "Sequence Name" ) +
                                         scale_x_continuous( breaks = c(seq( from = 1, to = Mosaic.Scale, by = Mosaic.Scale / 10 ), Mosaic.Scale), labels = c(bpX[seq( from = 1, to = Mosaic.Scale, by = Mosaic.Scale / 10 )], max(bpX)) ) + 
                                         scale_y_discrete( labels = c(ContigNames[3], ContigNames[2], ContigNames[1]))
                                       if( LabelBP == TRUE ){
                                         bars <- bars + theme( title = element_text(size = 14, colour = "black", face = "bold" ),
                                                               axis.text.y = element_text( colour=TickColour, size = TickSize ),
                                                               axis.title.y = element_text( size = LabelFontSize, colour = "black" ),
                                                               axis.text.x = element_text( colour=TickColour, size = TickSize),
                                                               axis.title.x = element_text( colour = "black", size = LabelFontSize)
                                         )
                                       } else {
                                         bars <- bars + theme( title = element_text(size = 14, colour = "black", face = "bold" ),
                                                               axis.text.y = element_text( colour=TickColour, size = TickSize ),
                                                               axis.title.y = element_text( size = LabelFontSize, colour = "black" ),
                                                               axis.text.x = element_blank(),
                                                               axis.title.x = element_blank()
                                         )
                                       }
                                       if( Title == T ) {
                                         bars <- bars + ggtitle( paste("Sequence Similarity Between Sequences for Triplet ", ContigNames[1], ":", ContigNames[2], ":", ContigNames[3], sep="") )
                                       }
                                       if( Legend == T ) {
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
                                         return( bars )
                                       } else {
                                         return( bars )
                                       }
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
                                     names(Thresholds) <- c( paste( ContigNames[1], ContigNames[2], sep=":" ), paste( ContigNames[1], ContigNames[2], sep=":" ), paste( ContigNames[2], ContigNames[3], sep=":" ) )
                                     cat( "Now beginning Block Search...\n\n" )
                                     Blocks <<- lapply( 1:3, function(i) block.find( SSTable[,c( 1:6, 6+i )], Thresholds[[i]] ) )
                                     names(Blocks) <<- names(Thresholds) <- c( paste( ContigNames[1], ContigNames[2], sep = ":" ), paste( ContigNames[1], ContigNames[3], sep=":" ), paste( ContigNames[2], ContigNames[3], sep=":" ) )
                                   },
                                 
                                 # Method for testing significance and dating of blocks.
                                 blockDate =
                                   function( dnaobj, mutation.rate = 10e-8, required.p = 0.005 ) {
                                     cat( "Now dating blocks" )
                                     ab.blocks <- lapply( Blocks[[1]], function(x) date.blocks( x, dnaobj, mutation.rate, 1, required.p ) )
                                     ac.blocks <- lapply( Blocks[[2]], function(x) date.blocks( x, dnaobj, mutation.rate, 2, required.p ) )
                                     bc.blocks <- lapply( Blocks[[3]], function(x) date.blocks( x, dnaobj, mutation.rate, 3, required.p ) )
                                     out.blocks <- list( ab.blocks, ac.blocks, bc.blocks )
                                     Blocks <<- mergeBandD( Blocks, out.blocks )
                                   },
                                 
                                 returnPair =
                                   function( sequence1, sequence2, data = T ) {
                                     pair <- c( which( ContigNames == sequence1 ), which( ContigNames == sequence2 ) )
                                     #pair <- c( which( c( SequenceA, SequenceB, SequenceC ) == sequence1 ), which( c( SequenceA, SequenceB, SequenceC ) == sequence2 ) )
                                     if( 1 %in% pair && 2 %in% pair ) {
                                       if( data == T ) {
                                         return( SSTable[,7] )
                                       } else {
                                         return( 1 )
                                       }
                                     } else {
                                       if( 1 %in% pair && 3 %in% pair ) {
                                         if( data ==T ) {
                                           return( SSTable[,8] )
                                         } else {
                                           return( 2 )
                                         }
                                       } else {
                                         if( 2 %in% pair && 3 %in% pair ) {
                                           if( data == T){
                                             return( SSTable[,9] )
                                           } else {
                                             return( 3 )
                                           }
                                         }
                                       }
                                     }
                                   }
                                 )
                               )