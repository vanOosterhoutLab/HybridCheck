# Plot the mosaic rainbow bars - Internal Function.
# Last altered by Ben J. Ward on 05/04/2013.

plot.mosaic <- function(x, xscale, condense, labfontsize, legfontsize) 
{
  # We can make the 3 colour references for the 3 different contigs:
  # Red - Contig A
  RefA <- expand.grid(contigb=seq(0, 100, by=1), contigc=seq(0, 100, by=1))
  RefA <- within(RefA, mix <- rgb(green=contigb, red=100, blue=contigc, maxColorValue=100))
  # Green - Contig B
  RefB <- expand.grid(contiga=seq(0, 100, by=1), contigc=seq(0, 100, by=1))
  RefB <- within(RefB, mix <- rgb(green=100, red=contiga, blue=contigc, maxColorValue=100))
  # Blue - Contig C
  RefC <- expand.grid(contiga=seq(0, 100, by=1), contigb=seq(0, 100, by=1))
  RefC <- within(RefC, mix <- rgb(green=contigb, red=contiga, blue=100, maxColorValue=100))
  
  # Generate some statistics to help with plotting data table formation.
  numberwin <- nrow(x[[1]]) # number of windows.
  yneeded <- (ceiling(numberwin/xscale)) # How many y tiles are needed to account for all windows.
  xAxis <- rep(c(1:xscale),yneeded) # Generate the xAxis variable.
  yAxis <- rep(1:yneeded, each = xscale) # Generate the yAxis variable.
  blanks <- (yneeded*xscale)-numberwin # There will be left overs, these are blanks, they only appear in the high res mosaic plot.
  
  # Create Matrices to contain all the values from which to work out the colours.
  ValuesA <- matrix(nrow=2, ncol=numberwin) # Values needed to work out the colours for red contig A.
  ValuesB <- matrix(nrow=2, ncol=numberwin) # Values needed to work out the colours for green contig B.
  ValuesC <- matrix(nrow=2, ncol=numberwin) # Values needed to work out the colours for blue contig C.
  
  # Fill out the three matrices from above using the distance data calculated.
  ValuesA[1,] <- as.numeric(x[[1]][,7])
  ValuesA[2,] <- as.numeric(x[[1]][,8])
  ValuesB[1,] <- as.numeric(x[[1]][,7])
  ValuesB[2,] <- as.numeric(x[[1]][,9])
  ValuesC[1,] <- as.numeric(x[[1]][,8]) 
  ValuesC[2,] <- as.numeric(x[[1]][,9])
  nopematrix <- matrix(nrow=2, ncol=blanks)
  
  # Append on the blank windows to the end of the colour mixing matrices...
  ValuesA <- cbind(ValuesA, nopematrix)
  ValuesB <- cbind(ValuesB, nopematrix)
  ValuesC <- cbind(ValuesC, nopematrix)
  # Make the plotting frames...
  Aframe <- data.frame(xaxis=xAxis, yaxis=yAxis, valuesB = ValuesA[1,], valuesC = ValuesA[2,])
  Bframe <- data.frame(xaxis=xAxis, yaxis=yAxis, valuesA = ValuesB[1,], valuesC = ValuesB[2,])
  Cframe <- data.frame(xaxis=xAxis, yaxis=yAxis, valuesA = ValuesC[1,], valuesB = ValuesC[2,])
  rm(ValuesA, ValuesB, ValuesC, nopematrix)
  
  if(condense==T){
    Aframe <- data.frame(xaxis=1:yneeded, yaxis=rep(1,yneeded), valuesB = as.vector(ceiling(tapply(Aframe$valuesB, Aframe$yaxis, mean))), valuesC = as.vector(ceiling(tapply(Aframe$valuesC, Aframe$yaxis, mean))))
    Bframe <- data.frame(xaxis=1:yneeded, yaxis=rep(1,yneeded), valuesA = as.vector(ceiling(tapply(Bframe$valuesA, Bframe$yaxis, mean))), valuesC = as.vector(ceiling(tapply(Bframe$valuesC, Bframe$yaxis, mean))))
    Cframe <- data.frame(xaxis=1:yneeded, yaxis=rep(1,yneeded), valuesA = as.vector(ceiling(tapply(Cframe$valuesA, Cframe$yaxis, mean))), valuesB = as.vector(ceiling(tapply(Cframe$valuesB, Cframe$yaxis, mean))))
    Aframe <- Aframe[complete.cases(Aframe),]
    Bframe <- Bframe[complete.cases(Bframe),]
    Cframe <- Cframe[complete.cases(Cframe),]
  }
  
  Aframe <- within(Aframe, mix <- apply(Aframe, 1, function(x) col.deter(as.numeric(x),RefA)))
  Bframe <- within(Bframe, mix <- apply(Bframe, 1, function(x) col.deter(as.numeric(x),RefB)))
  Cframe <- within(Cframe, mix <- apply(Cframe, 1, function(x) col.deter(as.numeric(x),RefC)))
  
  Aframe <- Aframe[,-c(3,4)]
  Bframe <- Bframe[,-c(3,4)]
  Cframe <- Cframe[,-c(3,4)]
  
  xwindowlabel<-paste("Sliding Window x",xscale,sep="")
  
  if(condense==T) {
    PlotA <- ggplot(Aframe, aes(x=xaxis, y=yaxis)) + 
      geom_raster(aes(fill=mix)) + 
      scale_fill_identity() +
      xlab(xwindowlabel) +
      ylab(paste("A).",x$ContigNames[1], sep=" ")) +
      theme(axis.text.y = element_blank(), axis.title.y = element_text(size=labfontsize, colour="black"), axis.ticks.y = element_blank(), legend.text = element_text(size=legfontsize))
    
    
    PlotB <- ggplot(Bframe, aes(x=xaxis, y=yaxis)) + 
      geom_raster(aes(fill=mix)) + 
      scale_fill_identity() +
      xlab(xwindowlabel) +
      ylab(paste("B).",x$ContigNames[2], sep=" ")) +
      theme(axis.text.y = element_blank(), axis.title.y = element_text(size=labfontsize, colour="black"), axis.ticks.y = element_blank(), legend.text = element_text(size=legfontsize))
    
    
    PlotC <- ggplot(Cframe, aes(x=xaxis, y=yaxis)) + 
      geom_raster(aes(fill=mix)) + 
      scale_fill_identity() +
      xlab(xwindowlabel) +
      ylab(paste("C).",x$ContigNames[3], sep=" ")) +
      theme(axis.text.y = element_blank(), axis.title.y = element_text(size=labfontsize, colour="black"), axis.ticks.y = element_blank(), legend.text = element_text(size=legfontsize))
    
  } else {
    
    PlotA <- ggplot(Aframe, aes(x=yaxis, y=xaxis)) + 
      geom_raster(aes(fill=mix)) + 
      scale_fill_identity() +
      ylab("Sliding Window (1:100)") +
      xlab("Sliding Window (x100)")
    
    PlotB <- ggplot(Bframe, aes(x=yaxis, y=xaxis)) + 
      geom_raster(aes(fill=mix)) + 
      scale_fill_identity() +
      ylab("Sliding Window (1:100)") +
      xlab("Sliding Window (x100)")
    
    PlotC <- ggplot(Cframe, aes(x=yaxis, y=xaxis)) + 
      geom_raster(aes(fill=mix)) + 
      scale_fill_identity() +
      ylab("Sliding Window (1:100)") +
      xlab(xwindowlabel)
  }
  
  return(as.HybRIDSbars(list(PlotA,PlotB,PlotC)))
  
}