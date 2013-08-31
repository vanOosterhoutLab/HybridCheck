plottingFrame <- data.frame( X = frame$X, Y = rep( c(3, 2, 1), each = mosaicscale), colour = c(frame$A_mix, frame$B_mix, frame$C_mix))
plot <- ggplot(plottingFrame, aes(x = X, y = as.factor(Y)) ) +
  geom_raster( aes( fill = colour, interpolate=T ) ) + scale_fill_identity() +
  xlab("Approximate Base Position") +
  ylab( "Sequence Name" ) +
  scale_x_continuous( breaks = c(seq( from = 1, to = mosaicscale, by = mosaicscale / 10 ), mosaicscale), labels = c(bpX[seq( from = 1, to = mosaicscale, by = mosaicscale / 10 )], max(bpX)) ) + 
  scale_y_discrete( labels = as.factor(c(Sequence3, Sequence2, Sequence1)) ) +
  theme( axis.text.y = element_text( colour="black", size = 12 ),
         axis.title.y = element_text( size = labfontsize, colour = "black" ),
         axis.text.x = element_text( colour="black", size = 12),
         axis.title.x = element_text( colour = "black", size = labfontsize)
         )