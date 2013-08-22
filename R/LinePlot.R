# Internal function for the generation of ggplot2 object of coloured lineplots of sequence similarity in the pairs of the DNA sequence triplets.
# Last Edited by Ben J. Ward on 30/04/2013.

plot.similarities <- function(x, leg, bpannotate, bpfreq, labfontsize, legfontsize, bpfixedxaxis) {
  combo <- unlist(lapply(combn(x$ContigNames, 2, simplify=FALSE), function(x) paste(x, collapse=":")))
  similarities <- as.matrix(x$Distances[,7:9])
  plotting.frame <- data.frame(basepos = rep(as.numeric(x$Distances[,4]), 3),
                               xrange = rep(c(1:nrow(similarities))),
                               yvalues = as.vector(similarities),
                               factors = rep(1:3, each = nrow(similarities)))
  if(bpannotate == TRUE && bpfixedxaxis == TRUE){
    linesplot <- ggplot(plotting.frame, aes(x=basepos, y=yvalues))
  } else {
    linesplot <- ggplot(plotting.frame, aes(x=xrange, y=yvalues))
  }
  linesplot <- linesplot +
    geom_line(aes(colour=factor(factors)), show_guide=leg, size=0.8) +
    ylim(0,100) +
    scale_colour_manual(name = "Pairwise Comparrisons", labels=c(combo[1], combo[2], combo[3]),values=c("yellow","purple","cyan")) +
    xlab("Sliding Window") + 
    ylab("% Sequence Similarity") +
    theme(axis.title.y = element_text(size=labfontsize, colour="black"), axis.title.x = element_text(size=labfontsize, colour="black"), axis.text.x = element_text(size=labfontsize, colour="black"), legend.text = element_text(size=legfontsize))
  if(bpannotate==TRUE && bpfixedxaxis == FALSE) {
    xlabsbreak <- c(1, seq(from = 0, to = nrow(x$Distances), by = bpfreq)[-1])
     linesplot <- linesplot +
       xlab("Sliding Window Central Position (actual bp)") +
       scale_x_continuous(breaks = c(xlabsbreak, nrow(x$Distances)), labels = c(as.numeric(plotting.frame$basepos[xlabsbreak]), as.numeric(plotting.frame$basepos[nrow(x$Distances)])))
  }
  return(as.HybRIDSlinesplot(linesplot))
}