# New plotting functions for absolute base scale

plot_mosaic <- function(ss, mosaic.scale, legfontsize = 12, labfontsize = 12) # ss is the SS object, mosaic.scale is how many vertical bars to put in the rainbow plot.
{
  # Now let's generate the reference colour palettes.
  RefA <- expand.grid(contigb = seq(0, 100, by = 1), contigc = seq(0, 100, by = 1))
  RefA <- within(RefA, mix <- rgb(green = contigb, red = 100, blue = contigc, maxColorValue = 100))
  RefB <- expand.grid(contiga = seq(0, 100, by = 1), contigc = seq(0, 100, by = 1))
  RefB <- within(RefB, mix <- rgb(green = 100, red = contiga, blue = contigc, maxColorValue = 100))
  RefC <- expand.grid(contiga = seq(0, 100, by = 1), contigb = seq(0, 100, by = 1))
  RefC <- within(RefC, mix <- rgb(green = contigb, red = contiga, blue = 100, maxColorValue = 100))
  
  # Now figure out the scale and data to go into each vertical bar:
  div <- ss$FullSeqLength / mosaic.scale
  sequence <- seq(from = 1, to = ss$FullSeqLength, by = div)
  seqend <- seq(from=div, to = ss$FullSeqLength, by = div)
  frame <- data.frame(bpstart = sequence, bpend = seqend)
  rm(sequence, seqend, div)
  
  # Now we go through each vertical bar, and for each one, we find the SSvalues to go in there, and we average them.
  AB <- round(apply(frame, 1, function(x) vertbar_create(ss, x, 7)))
  AC <- round(apply(frame, 1, function(x) vertbar_create(ss, x, 8)))
  BC <- round(apply(frame, 1, function(x) vertbar_create(ss, x, 9)))
  # We'll ust add them to the frame and remove the vectors.
  frame$AB <- AB
  frame$AC <- AC
  frame$BC <- BC
  frame$X <- 1:nrow(frame)
  rm(AB, AC, BC)
  
  frame$A_mix <- apply(frame, 1, function(x) col_deter( c( as.numeric(x[3]), as.numeric(x[4]) ), RefA))
  frame$B_mix <- apply(frame, 1, function(x) col_deter( c( as.numeric(x[3]), as.numeric(x[5]) ), RefB))
  frame$C_mix <- apply(frame, 1, function(x) col_deter( c( as.numeric(x[4]), as.numeric(x[5]) ), RefC))
  
  A_Bar <- ggplot(frame, aes(x = X, y = 1)) + 
    geom_raster(aes(fill = A_mix)) + scale_fill_identity() + ylab(paste("A).", ss$ContigNames[1], sep=" ")) +
    theme( axis.text.y = element_text(colour="white", size=4),
           axis.title.y = element_text( size = labfontsize, colour = "black"),
           axis.ticks.y = element_blank(),
           legend.text = element_text(size = legfontsize),
           axis.ticks.x = element_blank(),
           axis.title.x = element_blank(),
           axis.text.x = element_blank() )
  
  B_Bar <- ggplot(frame, aes(x = X, y = 1)) + 
    geom_raster(aes(fill = B_mix)) + scale_fill_identity() + ylab(paste("B).", ss$ContigNames[2], sep=" ")) +
    theme( axis.text.y = element_text(colour="white", size=4),
           axis.title.y = element_text( size = labfontsize, colour = "black"),
           axis.ticks.y = element_blank(),
           legend.text = element_text(size = legfontsize),
           axis.ticks.x = element_blank(),
           axis.title.x = element_blank(),
           axis.text.x = element_blank() )
  
  C_Bar <- ggplot(frame, aes(x = X, y = 1)) + 
    geom_raster(aes(fill = C_mix)) + scale_fill_identity() + ylab(paste("C).", ss$ContigNames[3], sep=" ")) +
    theme( axis.text.y = element_text(colour="white", size=4),
           axis.title.y = element_text( size = labfontsize, colour = "black"),
           axis.ticks.y = element_blank(),
           legend.text = element_text(size = legfontsize),
           axis.ticks.x = element_blank(),
           axis.title.x = element_blank(),
           axis.text.x = element_blank() )
  
  return(list(A_Bar, B_Bar, C_Bar))
}


# Internal function to create the vertical bars.
vertbar_create <- function(ssobj, plotframerow, whichcomp)
{
  bool1 <- as.numeric(ssobj$Distances[,5]) <= as.numeric(plotframerow[2])
  bool2 <- as.numeric(plotframerow[1]) <= as.numeric(ssobj$Distances[,6])
  index <- which(bool1 == bool2)
  return( mean(ssobj$Distances[index, whichcomp]) )
}


# Internal function to determine colours for the bars.
col_deter <- function(invalues, reference) 
{
  cols <- reference[reference[, 1] == as.numeric(invalues[1]) & reference[,2] == as.numeric(invalues[2]), 3]
  return(cols)
}