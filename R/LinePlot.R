# Internal function for the generation of ggplot2 object of coloured lineplots of sequence similarity in the pairs of the DNA sequence triplets.
# Last Edited by Ben J. Ward on 30/04/2013.

plot_similarities <- function(x, labfontsize, legfontsize) {
  combo <- unlist(lapply(combn(x$ContigNames, 2, simplify=FALSE), function(x) paste(x, collapse=":")))
  similarities <- as.matrix(x$Distances[,7:9])
  plotting.frame <- data.frame(basepos = rep(as.numeric(x$Distances[,4]), 3),
                               xrange = rep(c(1:nrow(similarities))),
                               yvalues = as.vector(similarities),
                               factors = rep(1:3, each = nrow(similarities)))
  
  linesplot <- ggplot(plotting.frame, aes(x=basepos, y=yvalues)) + geom_line(aes(colour=factor(factors)), show_guide=T, size=0.8) +
    ylim(0,100) + 
    scale_colour_manual(name = "Pairwise Comparrisons", labels=c(combo[1], combo[2], combo[3]),values=c("yellow","purple","cyan")) +
    xlab("Base Position") +
    ylab("% Sequence Similarity") +
    theme(axis.title.y = element_text(size=labfontsize, colour="black"),
          axis.title.x = element_text(size=labfontsize, colour="black"), 
          axis.text.x = element_text(size=labfontsize, colour="black"), 
          legend.text = element_text(size=legfontsize) )
    
  return(linesplot)
}