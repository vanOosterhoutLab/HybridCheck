# An internal function tha applies all the plotting parameters to a plot - mostly to save code typing!
applyPlottingParams <- function( plot, parameters, title = "") {
  if( parameters$PlotTitle == TRUE ) {
    plot <- plot + 
      ggtitle( title ) +
      theme( title = element_text(size = parameters$TitleSize, colour = parameters$TitleColour, face = parameters$TitleFace ),
             legend.text = element_text( size = parameters$LegendFontSize ) )
  } else {
    plot <- plot + theme( title = element_blank(),
                          legend.text = element_text( size = parameters$LegendFontSize ) )
  }
  if( parameters$XTitle == TRUE ) {
    plot <- plot + theme( axis.title.x = element_text( size = parameters$XTitleSize, colour = parameters$XTitleColour ) )
  } else {
    plot <- plot + theme( axis.title.x = element_blank() )
  }
  if( parameters$YTitle == TRUE ) {
    plot <- plot + theme( axis.title.y = element_text( size = parameters$YTitleFontSize, colour = parameters$YTitleColour ) )
  } else {
    plot <- plot + theme( axis.title.y = element_blank() )
  }
  if( parameters$XLabel == TRUE ) {
    plot <- plot + theme( axis.text.x = element_text( size = parameters$XLabelSize, colour = parameters$XLabelColour ) )
  } else {
    plot <- plot + theme( axis.text.x = element_blank() )
  }
  if( parameters$YLabel == TRUE ) {
    plot <- plot + theme( axis.text.y = element_text( size = parameters$YLabelSize, colour = parameters$YLabelColour ) )
  } else {
    plot <- plot + theme( axis.text.y = element_blank() )
  }
  return(plot)
}