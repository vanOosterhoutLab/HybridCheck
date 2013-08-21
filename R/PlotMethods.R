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
col.deter <- function(invalues, reference) {
  cols <- reference[reference[,1]==invalues[2] & reference[,2]==invalues[3],3]
  return(cols)
}


# Internal function for subsetting the HybRIDS data objects.
subseq <- function(x, s1, s2){
  contains <- unlist(lapply(x, function(n) s1 %in% n$ContigNames && s2 %in% n$ContigNames))
  if(!TRUE %in% contains){
    stop("Those contig names are not in the present dataset.")
  }
  subx <- x[contains]
  return(subx)
}


# Internal function for getting the right SS values.
subsim <- function(input, p){
  if(all(p == c(1,2))){
    return(input$Distance[,7])
  } else {
    if(all(p == c(2,3))){
      return(input$Distance[,8])
    } else {
      if(all(p == c(1,3))){
        return(input$Distance[,9])
      }
    }
  }
}


# # Internal function for making the data fit an absolute bp by bp scale.
# absoluteBP <- function(x){
#   bpdf <- data.frame(matrix(ncol=2, nrow = subsetx[[1]]$FullSeqLength ))
#   windowcover <- lapply(1:nrow(x$Distance), function(i) c(x$Distance$ActualStart[i]:x$Distance$ActualEnd[i]))
#   
#   
#   
#   
#   
#   
# }





# Plotting method for a set of triplet analyses.
plot.HybRIDSseqsimSET <- function(x, Sequence1, Sequence2) {
  subsetx <- subseq(x, Sequence1, Sequence2)
  cat("Check1")
  pair <- lapply(subsetx, function(x) c(which(x$ContigNames==Sequence1), which(x$ContigNames==Sequence2)))
  cat("Check2")
  dflength <- sum(unlist(lapply(subsetx, function(x) nrow(x$Distances))))
  cat("Check3")
  plotting.frame <- data.frame(matrix(nrow = dflength, ncol = 9))
  names(plotting.frame) <- c("WindowCenter", "WindowStart", "WindowEnd", "ActualCenter", "ActualStart", "ActualEnd", "SSVals", "TripletSet", "xvals")
  plotting.frame$xvals <- unlist(lapply(subsetx, function(x) 1:nrow(x$Distance)))
  plotting.frame$WindowCenter <- unlist(lapply(subsetx, function(x) x$Distance$WindowCenter))
  plotting.frame$WindowStart <- unlist(lapply(subsetx, function(x) x$Distance$WindowStart))
  plotting.frame$WindowEnd <- unlist(lapply(subsetx, function(x) x$Distance$WindowEnd))
  plotting.frame$ActualCenter <- unlist(lapply(subsetx, function(x) x$Distance$ActualCenter))
  plotting.frame$ActualStart <- unlist(lapply(subsetx, function(x) x$Distance$ActualStart))
  plotting.frame$ActualEnd <- unlist(lapply(subsetx, function(x) x$Distance$ActualEnd))
  plotting.frame$SSVals <- unlist(lapply(1:length(subsetx), function(i) subsim(subsetx[[i]],pair[[i]])))
  plotting.frame$TripletSet <- as.factor(unlist(lapply(subsetx, function(x) rep(paste(x$ContigNames, collapse=":"), nrow(x$Distances)))))
  linesplot <- ggplot(plotting.frame, aes(x=ActualCenter, y=SSVals)) +
    geom_line(aes(colour=TripletSet), show_guide=T, size=0.8) +
    ylab("% Sequence Similarity") +
    xlab("Base Position")
  return(linesplot)
}


dfdates <- function(){
  
}


# Plotting method for the location of blocks and the date of the blocks.
# plot.HybRIDSblocks <- function(input) {
#   plotting.frame <- data.frame(matrix(nrow = sum(unlist(lapply(input[[2]], function(x) unlist(lapply(x, function(y) nrow(y)))))), ncol=7))
#   names(plotting.frame) <- c("Pair", "Threshold", "Five", "Fifty", "Size", "SNPs", "SSVals", "TripletSet", "xvals")
#   
#   
# }
# 
# 
# plot.HybRIDSdatedblocks <- function(){
#   
# }

