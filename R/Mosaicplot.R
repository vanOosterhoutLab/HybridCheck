# Plot the mosaic rainbow bars - Internal Function.

plot.mosaic <- function (x, xscale, condense, labfontsize, legfontsize) 
{
  RefA <- expand.grid(contigb = seq(0, 100, by = 1), contigc = seq(0, 
                                                                   100, by = 1))
  RefA <- within(RefA, mix <- rgb(green = contigb, red = 100, 
                                  blue = contigc, maxColorValue = 100))
  RefB <- expand.grid(contiga = seq(0, 100, by = 1), contigc = seq(0, 
                                                                   100, by = 1))
  RefB <- within(RefB, mix <- rgb(green = 100, red = contiga, 
                                  blue = contigc, maxColorValue = 100))
  RefC <- expand.grid(contiga = seq(0, 100, by = 1), contigb = seq(0, 
                                                                   100, by = 1))
  RefC <- within(RefC, mix <- rgb(green = contigb, red = contiga, 
                                  blue = 100, maxColorValue = 100))
  numberwin <- nrow(x[[1]])
  yneeded <- (ceiling(numberwin/xscale))
  xAxis <- rep(c(1:xscale), yneeded)
  yAxis <- rep(1:yneeded, each = xscale)
  blanks <- (yneeded * xscale) - numberwin
  ValuesA <- matrix(nrow = 2, ncol = numberwin)
  ValuesB <- matrix(nrow = 2, ncol = numberwin)
  ValuesC <- matrix(nrow = 2, ncol = numberwin)
  ValuesA[1, ] <- as.numeric(x[[1]][, 7])
  ValuesA[2, ] <- as.numeric(x[[1]][, 8])
  ValuesB[1, ] <- as.numeric(x[[1]][, 7])
  ValuesB[2, ] <- as.numeric(x[[1]][, 9])
  ValuesC[1, ] <- as.numeric(x[[1]][, 8])
  ValuesC[2, ] <- as.numeric(x[[1]][, 9])
  nopematrix <- matrix(nrow = 2, ncol = blanks)
  ValuesA <- cbind(ValuesA, nopematrix)
  ValuesB <- cbind(ValuesB, nopematrix)
  ValuesC <- cbind(ValuesC, nopematrix)
  Aframe <- data.frame(xaxis = xAxis, yaxis = yAxis, valuesB = ValuesA[1, 
                                                                       ], valuesC = ValuesA[2, ])
  Bframe <- data.frame(xaxis = xAxis, yaxis = yAxis, valuesA = ValuesB[1, 
                                                                       ], valuesC = ValuesB[2, ])
  Cframe <- data.frame(xaxis = xAxis, yaxis = yAxis, valuesA = ValuesC[1, 
                                                                       ], valuesB = ValuesC[2, ])
  rm(ValuesA, ValuesB, ValuesC, nopematrix)
  if (condense == T) {
    Aframe <- data.frame(xaxis = 1:yneeded, yaxis = rep(1, 
                                                        yneeded), valuesB = as.vector(ceiling(tapply(Aframe$valuesB, 
                                                                                                     Aframe$yaxis, mean))), valuesC = as.vector(ceiling(tapply(Aframe$valuesC, 
                                                                                                                                                               Aframe$yaxis, mean))))
    Bframe <- data.frame(xaxis = 1:yneeded, yaxis = rep(1, 
                                                        yneeded), valuesA = as.vector(ceiling(tapply(Bframe$valuesA, 
                                                                                                     Bframe$yaxis, mean))), valuesC = as.vector(ceiling(tapply(Bframe$valuesC, 
                                                                                                                                                               Bframe$yaxis, mean))))
    Cframe <- data.frame(xaxis = 1:yneeded, yaxis = rep(1, 
                                                        yneeded), valuesA = as.vector(ceiling(tapply(Cframe$valuesA, 
                                                                                                     Cframe$yaxis, mean))), valuesB = as.vector(ceiling(tapply(Cframe$valuesB, 
                                                                                                                                                               Cframe$yaxis, mean))))
    Aframe <- Aframe[complete.cases(Aframe), ]
    Bframe <- Bframe[complete.cases(Bframe), ]
    Cframe <- Cframe[complete.cases(Cframe), ]
  }
  Aframe <- within(Aframe, mix <- apply(Aframe, 1, function(x) col.deter(as.numeric(x), 
                                                                         RefA)))
  Bframe <- within(Bframe, mix <- apply(Bframe, 1, function(x) col.deter(as.numeric(x), 
                                                                         RefB)))
  Cframe <- within(Cframe, mix <- apply(Cframe, 1, function(x) col.deter(as.numeric(x), 
                                                                         RefC)))
  Aframe <- Aframe[, -c(3, 4)]
  Bframe <- Bframe[, -c(3, 4)]
  Cframe <- Cframe[, -c(3, 4)]
  xwindowlabel <- paste("Sliding Window x", xscale, sep = "")
  if (condense == T) {
    PlotA <- ggplot(Aframe, aes(x = xaxis, y = yaxis)) + 
      geom_raster(aes(fill = mix)) + scale_fill_identity() + 
      xlab(xwindowlabel) + ylab(paste("A).", x$ContigNames[1], 
                                      sep = " ")) + theme(axis.text.y = element_blank(), 
                                                          axis.title.y = element_text(size = labfontsize, colour = "black"), 
                                                          axis.ticks.y = element_blank(), legend.text = element_text(size = legfontsize))
    PlotB <- ggplot(Bframe, aes(x = xaxis, y = yaxis)) + 
      geom_raster(aes(fill = mix)) + scale_fill_identity() + 
      xlab(xwindowlabel) + ylab(paste("B).", x$ContigNames[2], 
                                      sep = " ")) + theme(axis.text.y = element_blank(), 
                                                          axis.title.y = element_text(size = labfontsize, colour = "black"), 
                                                          axis.ticks.y = element_blank(), legend.text = element_text(size = legfontsize))
    PlotC <- ggplot(Cframe, aes(x = xaxis, y = yaxis)) + 
      geom_raster(aes(fill = mix)) + scale_fill_identity() + 
      xlab(xwindowlabel) + ylab(paste("C).", x$ContigNames[3], 
                                      sep = " ")) + theme(axis.text.y = element_blank(), 
                                                          axis.title.y = element_text(size = labfontsize, colour = "black"), 
                                                          axis.ticks.y = element_blank(), legend.text = element_text(size = legfontsize))
  }
  else {
    PlotA <- ggplot(Aframe, aes(x = yaxis, y = xaxis)) + 
      geom_tile(aes(fill = mix), color = "white") + scale_fill_identity() + 
      ylab("Sliding Window (1:100)") + xlab("Sliding Window (x100)")
    PlotB <- ggplot(Bframe, aes(x = yaxis, y = xaxis)) + 
      geom_tile(aes(fill = mix), color = "white") + scale_fill_identity() + 
      ylab("Sliding Window (1:100)") + xlab("Sliding Window (x100)")
    PlotC <- ggplot(Cframe, aes(x = yaxis, y = xaxis)) + 
      geom_tile(aes(fill = mix), color = "white") + scale_fill_identity() + 
      ylab("Sliding Window (1:100)") + xlab(xwindowlabel)
  }
  return(as.HybRIDSbars(list(PlotA, PlotB, PlotC)))
}