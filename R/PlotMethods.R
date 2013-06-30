# Methods for plotting the various HybRIDS classes.
# Created by Ben J. Ward 30/04/2013.
# Last Edited by Ben J. Ward on 30/04/2013.

#' Method for Plotting objects of the class HybRIDSseqsim, associated with the HybRIDS package.
#' 
#' @S3method plot HybRIDSseqsim
#' @method plot HybRIDSseqsim
plot.HybRIDSseqsim <- function(x, linesplot=T, legends=T, baseannotate=T, bpfreq=500, densityplot=F, mosaic.bars=T, mosaic.scale=50, condense.mosaics=T, labfontsize=14, legfontsize=14, onecanvas=T) {
  if(linesplot==F && mosaic.bars==F && densityplot==F) stop("You have specified no type of plot to be plotted.
                                          You need to set one of these to true.")
  if(linesplot==T){
    LinesPlot <- plot.similarities(x, legends, baseannotate, bpfreq, labfontsize,legfontsize)
    if(onecanvas==F) plot(LinesPlot)
  }
  if(mosaic.bars==T){
    BarsPlot <- plot.mosaic(x, mosaic.scale, condense.mosaics, labfontsize,legfontsize)
    if(onecanvas==F){
      if(legends==T){
        print(consolidate.bars(BarsPlot))
      } else {
        print(consolidate.bars(BarsPlot,pyramid=F))
      }
    } 
  }
  if(densityplot==T){
    DensityFrame <- data.frame(comp=rep(1:3, each=nrow(x$Distances)), SS=c(x$Distances[,7],x$Distances[,8],x$Distances[,9]), SS.mean=rep(c(mean(x$Distances[,7]),mean(x$Distances[,8]),mean(x$Distances[,9])), each=nrow(x$Distances)), half.SD=rep(c(mean(x$Distances[,7])+(sd(x$Distances[,7])/2),mean(x$Distances[,8])+(sd(x$Distances[,8])/2),mean(x$Distances[,9])+(sd(x$Distances[,9])/2)), each=nrow(x$Distances)))
    DensityPlot <- ggplot(DensityFrame, aes(x=SS, colour=factor(comp))) + 
      geom_density() + 
      geom_vline(data=DensityFrame, aes(xintercept=SS.mean,  colour=factor(comp)), linetype="solid", size=1) +
      geom_vline(data=DensityFrame, aes(xintercept=half.SD,  colour=factor(comp)), linetype="dashed", size=1) +
      scale_colour_manual(name = "Pairwise Comparrisons", labels=c(paste(x$ContigNames[1],x$ContigNames[2],sep=":"),paste(x$ContigNames[1],x$ContigNames[3],sep=":"),paste(x$ContigNames[2],x$ContigNames[3],sep=":")),values=c("yellow","purple","cyan"))
    if(onecanvas==F){
      if(legends==T){
        plot(DensityPlot)
      } else {
        plot(DensityPlot + theme(legend.position="none", axis.text.y=element_text(size=10), axis.ticks.y=element_blank()))
      }
    } 
  }
  
  if(onecanvas==T && length(which(c(exists("BarsPlot"),exists("LinesPlot"),exists("DensityPlot"))==TRUE))==1){
    stop("You have selected to plot only one plot, but to combine multiple plots one one canvas - your command conflicts.\n either plot more than one plot or turn the onecanvas option to false.")
  }
  
  if(onecanvas==T){
    if(exists("BarsPlot")==F && exists("LinesPlot")==F && exists("DensityPlot")==F){
      stop("ERROR: You have plotted no plots!")
    } else {
      if(exists("BarsPlot")==T && exists("LinesPlot")==T && exists("DensityPlot")==F){
        BarsPlot[[1]] <- BarsPlot[[1]]+
          scale_y_continuous(breaks=1, labels="100") +
          theme(axis.text.y = element_text(size=10, colour="white"), axis.title.y = element_text(), axis.ticks.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
        
        BarsPlot[[2]] <- BarsPlot[[2]]+
          scale_y_continuous(breaks=1, labels= "100") +
          theme(axis.text.y = element_text(size=10, colour="white"), axis.title.y = element_text(), axis.ticks.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
        
        BarsPlot[[3]] <- BarsPlot[[3]]+
          scale_y_continuous(breaks=1, labels="100") +
          theme(axis.text.y = element_text(size=10, colour="white"), axis.title.y = element_text(), axis.ticks.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
        
        tmp <- ggplot_gtable(ggplot_build(LinesPlot))
        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        if(length(leg)>0){
          legend <- tmp$grobs[[leg]]
          grid.arrange( arrangeGrob(consolidate.bars(BarsPlot),
                                    arrangeGrob(LinesPlot + theme(legend.position="none", axis.text.y=element_text(size=10), axis.ticks.y=element_blank()),
                                                legend, widths=c(1,0.15), ncol=2), ncol=1))
        } else {
          grid.arrange( arrangeGrob(consolidate.bars(BarsPlot,pyramid=F), widths=c(1,0.15),
                                    LinesPlot + theme(legend.position="none", axis.text.y=element_text(size=10), axis.ticks.y=element_blank()), ncol=1))
        }
        
      } else {
        
        # When the user wants Bars, Lines & The density charts one one canvas...
        
        if(exists("BarsPlot")==T && exists("LinesPlot")==T && exists("DensityPlot")==T){
          BarsPlot[[1]] <- BarsPlot[[1]]+
            scale_y_continuous(breaks=1, labels="100") +
            theme(axis.text.y = element_text(size=10, colour="white"), axis.title.y = element_text(), axis.ticks.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
          
          BarsPlot[[2]] <- BarsPlot[[2]]+
            scale_y_continuous(breaks=1, labels= "100") +
            theme(axis.text.y = element_text(size=10, colour="white"), axis.title.y = element_text(), axis.ticks.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
          
          BarsPlot[[3]] <- BarsPlot[[3]]+
            scale_y_continuous(breaks=1, labels="100") +
            theme(axis.text.y = element_text(size=10, colour="white"), axis.title.y = element_text(), axis.ticks.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
          
          tmp1 <- ggplot_gtable(ggplot_build(LinesPlot))
          leg1 <- which(sapply(tmp1$grobs, function(x) x$name) == "guide-box")
          
          if(length(leg1)>0){
            legend1 <- tmp1$grobs[[leg1]]
            tmp2 <- ggplot_gtable(ggplot_build(DensityPlot))
            leg2 <- which(sapply(tmp2$grobs, function(x) x$name) == "guide-box")
            legend2 <- tmp2$grobs[[leg2]]
            grid.arrange( arrangeGrob(consolidate.bars(BarsPlot),
                                      arrangeGrob(LinesPlot + theme(legend.position="none", axis.text.y=element_text(size=10), axis.ticks.y=element_blank()),
                                                  legend1, widths=c(1,0.15), ncol=2), arrangeGrob(DensityPlot + theme(legend.position="none", axis.text.y=element_text(size=10), axis.ticks.y=element_blank()),
                                                                                                  legend2, widths=c(1,0.15), ncol=2), ncol=1))
          } else {
            grid.arrange( arrangeGrob(consolidate.bars(BarsPlot,pyramid=F), widths=c(1,0.15),
                                      arrangeGrob(LinesPlot + theme(legend.position="none", axis.text.y=element_text(size=10), axis.ticks.y=element_blank()), DensityPlot + theme(legend.position="none", axis.text.y=element_text(size=10), axis.ticks.y=element_blank()), ncol=1)))
          }
        } else {
          if(exists("BarsPlot")==T && exists("LinesPlot")==F && exists("DensityPlot")==T){
            BarsPlot[[1]] <- BarsPlot[[1]]+
              scale_y_continuous(breaks=1, labels="100") +
              theme(axis.text.y = element_text(size=10, colour="white"), axis.title.y = element_text(), axis.ticks.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
            
            BarsPlot[[2]] <- BarsPlot[[2]]+
              scale_y_continuous(breaks=1, labels= "100") +
              theme(axis.text.y = element_text(size=10, colour="white"), axis.title.y = element_text(), axis.ticks.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
            
            BarsPlot[[3]] <- BarsPlot[[3]]+
              scale_y_continuous(breaks=1, labels="100") +
              theme(axis.text.y = element_text(size=10, colour="white"), axis.title.y = element_text(), axis.ticks.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
            
            
            if(legends==T){
              tmp <- ggplot_gtable(ggplot_build(DensityPlot))
              leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
              legend <- tmp$grobs[[leg]]
              grid.arrange( arrangeGrob(consolidate.bars(BarsPlot),
                                        arrangeGrob(DensityPlot + theme(legend.position="none", axis.text.y=element_text(size=10), axis.ticks.y=element_blank()),
                                                    legend, widths=c(1,0.15), ncol=2), ncol=1))
            } else {
              grid.arrange( consolidate.bars(BarsPlot,pyramid=F),
                            DensityPlot + theme(legend.position="none", axis.text.y=element_text(size=10), axis.ticks.y=element_blank()), ncol=1)
            }
          } else {
            if(exists("BarsPlot")==F && exists("LinesPlot")==T && exists("DensityPlot")==T){
              tmp1 <- ggplot_gtable(ggplot_build(LinesPlot))
              leg1 <- which(sapply(tmp1$grobs, function(x) x$name) == "guide-box")
              
              if(length(leg1)>0){
                legend1 <- tmp1$grobs[[leg1]]
                tmp2 <- ggplot_gtable(ggplot_build(DensityPlot))
                leg2 <- which(sapply(tmp2$grobs, function(x) x$name) == "guide-box")
                legend2 <- tmp2$grobs[[leg2]]
                grid.arrange( arrangeGrob(arrangeGrob(LinesPlot + theme(legend.position="none", axis.text.y=element_text(size=10), axis.ticks.y=element_blank()),
                                                      legend1, widths=c(1,0.15), ncol=2), arrangeGrob(DensityPlot + theme(legend.position="none", axis.text.y=element_text(size=10), axis.ticks.y=element_blank()),
                                                                                                      legend2, widths=c(1,0.15), ncol=2), ncol=1))
              } else {
                grid.arrange(LinesPlot + theme(legend.position="none", axis.text.y=element_text(size=10), axis.ticks.y=element_blank()), DensityPlot + theme(legend.position="none", axis.text.y=element_text(size=10), axis.ticks.y=element_blank()), ncol=1)
              } 
            }
          }
        }
      }
    }
  }
}
      



consolidate.bars <- function(x, pyramid=T){
  legend <- readPNG(system.file("extdata/rgblegend.png", package="HybRIDS"), TRUE)
  if (names(dev.cur()) == "windows") {
    # windows device doesn’t support semi-transparency so we’ll need
    # to flatten the image
    transparent <- legend[,,4] == 0
    legend <- as.raster(legend[,,1:3])
    legend[transparent] <- NA
  }
  legendgrob <- rasterGrob(image=legend)
  if(pyramid==T){
    barsgrob <- arrangeGrob(arrangeGrob(x[[1]],
                            x[[2]],
                            x[[3]],ncol=1), legendgrob, widths=c(1,0.15), ncol=2)
    return(barsgrob)
  } else {
    if(pyramid==F){
      barsgrob <- arrangeGrob(x[[1]],
                              x[[2]],
                              x[[3]],ncol=1)
      return(barsgrob)
    } 
  }
}


print.HybRIDSbars <- function(x){
  y <- plot(x)
  return(y)
}



# Function to Determine the colour to use for the mosaic rainbow plots.
col.deter <- function(invalues,reference) {
  if(is.na(invalues[3]) && is.na(invalues[4])){
    cols <- "#000000"
  } else {
    cols <- reference[reference[,1]==invalues[3] & reference[,2]==invalues[4],3]
  }
  return(cols)
}



# plot.HybRIDSseqsimSET <- function(x) {
#   plottingframe <- matrix(nrow=nrow(x[[1]]$Distances), ncol=length(x)*3)
#   for(i in 1:length(x)){
#     indicies <- c(1,2,3)+(3*(i-1))
#     plottingframe[,indicies] <- x[[i]]$Distances[,c(7,8,9)]
#     colnames(plottingframe)[,indicies] <- "yo"
#   }
#   
#   
#   
#   
# }