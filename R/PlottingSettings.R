PlottingSettings <- setRefClass("PlottingSettings",
                                
                                fields = list(
                                  What = "character",
                                  PlotTitle = "logical",
                                  CombinedTitle = "logical",
                                  TitleSize = "numeric",
                                  TitleFace = "character",
                                  TitleColour = "character",
                                  XLabels = "logical",
                                  YLabels = "logical",
                                  XTitle = "logical",
                                  XTitleFontSize = "numeric",
                                  XTitleColour = "character",
                                  XLabelSize = "numeric",
                                  XLabelColour = "character",
                                  YTitle = "logical",
                                  YTitleFontSize = "numeric",
                                  YTitleColour = "character",
                                  YLabelSize = "numeric",
                                  YLabelColour = "character",
                                  Legends = "logical",
                                  LegendFontSize = "numeric",
                                  MosaicScale = "numeric"),
                                
                                methods = list(
                                  initialize = function(){
                                    What <<- c("Bars", "Lines")
                                    PlotTitle <<- TRUE
                                    CombinedTitle <<- FALSE 
                                    TitleSize <<- 14
                                    TitleFace <<- "bold"
                                    TitleColour <<- "black"
                                    XLabels <<- TRUE
                                    YLabels <<- TRUE
                                    XTitle <<- TRUE
                                    XTitleFontSize <<- 12
                                    XTitleColour <<- "black"
                                    XLabelSize <<- 10
                                    XLabelColour <<- "black"
                                    YTitle <<- TRUE
                                    YTitleFontSize <<- 12
                                    YTitleColour <<- "black"
                                    YLabelSize <<- 10
                                    YLabelColour <<- "black"
                                    Legends <<- TRUE
                                    LegendFontSize <<- 12
                                    MosaicScale <<- 500
                                  },
                                  
                                  setWhat = function(value){
                                    if(!all(value == "Lines" | value == "Bars") || length(value) > 2){stop("Input value must be 'Bars' and/or 'Lines'")}
                                    What <<- value
                                  },
                                  
                                  textSummary =
                                    function(){
                                      return(paste('Settings for plotting recombination signal in triplets:\n',
                                                   '-------------------------------------------------------\n\n',
                                                   'What graphs will be plotted (What): ', paste(What, collapse=", "),
                                                   '\n\nWhether to include a title with the plot. (PlotTitle): ',
                                                   PlotTitle,
                                                   '\n\nOne overall title for a Lines and Bars plot on the same canvas (CombinedTitle): ',
                                                   CombinedTitle,
                                                   '\n\nFont size of titles (TitleSize): ', TitleSize,
                                                   '\n\nFace of the title (TitleFace): ', TitleFace,
                                                   '\n\nColour of the plot titles (TitleColour): ', TitleColour,
                                                   '\n\nInclude labels of the x axes in plots (XLabels): ', XLabels,
                                                   '\n\nSize of the x axes labels (XLabelSize): ', XLabelSize,
                                                   '\n\nColiur of the x axes labels (XLabelColour): ', XLabelColour,
                                                   '\n\nInclude the title of the x axis in plots (XTitle): ', XTitle,
                                                   '\n\nFont size of the x axis title (XTitleFontSize): ', XTitleFontSize,
                                                   '\n\nColour of the x axis title (XTitleColour): ', XTitleColour,
                                                   '\n\nInclude labels of the y axes in plots (YLabels): ', YLabels,
                                                   '\n\nSize of the y axes labels (YLabelSize): ', YLabelSize,
                                                   '\n\nColiur of the y axes labels (YLabelColour): ', YLabelColour,
                                                   '\n\nInclude the title of the y axis in plots (YTitle): ', YTitle,
                                                   '\n\nFont size of the y axis title (YTitleFontSize): ', YTitleFontSize,
                                                   '\n\nColour of the y axis title (YTitleColour): ', YTitleColour,
                                                   '\n\nInclude legends for plots (Legends): ', Legends,
                                                   '\n\nFont size for legends (LegendFontSize): ', LegendFontSize,
                                                   '\n\nMosaicScale: ', MosaicScale      
                                        ))
                                    },
                                  
                                  show = function(){
                                    cat(textSummary())
                                  },
                                  
                                  setSettings = function(...){
                                    settings <- list(...)
                                    parameters <- names(settings)
                                    for(i in 1:length(settings)){
                                      switch(
                                        parameters[i],
                                        What = setWhat(settings[[i]]),
                                        PlotTitle = PlotTitle <<- settings[[i]],
                                        CombinedTitle = CombinedTitle <<- settings[[i]],
                                        TitleSize = TitleSize <<- settings[[i]],
                                        TitleFace = TitleFace <<- settings[[i]],
                                        TitleColour = TitleColour <<- settings[[i]],
                                        XLabels = XLabels <<- settings[[i]],
                                        YLabels = YLabels <<- settings[[i]],
                                        XTitle = XTitle <<- settings[[i]],
                                        XTitleFontSize = XTitleFontSize <<- settings[[i]],
                                        XTitleColour = XTitleColour <<- settings[[i]],
                                        XLabelSize = XLabelSize <<- settings[[i]],
                                        XLabelColour = XLabelColour <<- settings[[i]],
                                        YTitle = YTitle <<- settings[[i]],
                                        YTitleFontSize = YTitleFontSize <<- settings[[i]],
                                        YTitleColour = YTitleColour <<- settings[[i]],
                                        YLabelSize = YLabelSize <<- settings[[i]],
                                        YLabelColour = YLabelColour <<- settings[[i]],
                                        Legends = Legends <<- settings[[i]],
                                        LegendFontSize = LegendFontSize <<- settings[[i]],
                                        MosaicScale = MosaicScale <<- settings[[i]]
                                        )
                                    }
                                  }
                                  )
                                )

applyPlottingParams <- function(plot, parameters, title = ""){
  if(parameters$PlotTitle == TRUE){
    plot <- plot + 
      ggtitle(title) +
      theme(title = element_text(size = parameters$TitleSize, colour = parameters$TitleColour, face = parameters$TitleFace),
             legend.text = element_text(size = parameters$LegendFontSize))
  } else {
    plot <- plot + theme(title = element_blank(),
                          legend.text = element_text(size = parameters$LegendFontSize))
  }
  if(parameters$XTitle == TRUE){
    plot <- plot + theme(axis.title.x = element_text(size = parameters$XTitleFontSize, colour = parameters$XTitleColour))
  } else {
    plot <- plot + theme(axis.title.x = element_blank())
  }
  if(parameters$YTitle == TRUE){
    plot <- plot + theme(axis.title.y = element_text(size = parameters$YTitleFontSize, colour = parameters$YTitleColour))
  } else {
    plot <- plot + theme(axis.title.y = element_blank())
  }
  if(parameters$XLabels == TRUE){
    plot <- plot + theme(axis.text.x = element_text(size = parameters$XLabelSize, colour = parameters$XLabelColour))
  } else {
    plot <- plot + theme(axis.text.x = element_blank())
  }
  if(parameters$YLabels == TRUE){
    plot <- plot + theme(axis.text.y = element_text(size = parameters$YLabelSize, colour = parameters$YLabelColour))
  } else {
    plot <- plot + theme(axis.text.y = element_blank())
  }
  return(plot)
}