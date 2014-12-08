library(grid)

rowMax.units <- function(u, nrow){ # rowMax with a fake matrix of units
  matrix.indices <- matrix(seq_along(u), nrow=nrow)
  do.call(unit.c, lapply(seq(1, nrow), function(ii) {
   max(u[matrix.indices[ii, ]])
  }))
}

colMax.units <- function(u, ncol){ # colMax with a fake matrix of units
  matrix.indices <- matrix(seq_along(u), ncol=ncol)
  do.call(unit.c, lapply(seq(1, ncol), function(ii) {
   max(u[matrix.indices[, ii]])
  }))
}


makeTableGrobs <- function(e, ncol, nrow, equal.width = F, equal.height=F, 
        just = c("center", "center"),
        padding.h = unit(4, "mm"), padding.v=unit(4, "mm"), 
        gpar.text = gpar(col="black", cex=1),
        gpar.fill = gpar(fill = "grey95", col="white")) {

 n <- length(e) # number of labels

 if(missing(ncol) & missing(nrow)){
 nm <- n2mfrow(n)      # pretty default layout
 ncol = nm[1]
 nrow = nm[2]
}

makeOneLabel <- function(label.ind){
textGrob(label=e[label.ind], gp=gpar.text, name=paste("cells-label-",label.ind, sep=""))
}

makeOneCell <- function(label.ind){
rectGrob(gp=gpar.fill, name=paste("cells-fill-",label.ind, sep=""))
}


  lg <- lapply(seq_along(e), makeOneLabel) # list of text grobs
  lf <- lapply(seq_along(e), makeOneCell) # list of rect grobs
  
  wg <- lapply(lg, grobWidth) # list of grob widths
  hg <- lapply(lg, grobHeight) # list of grob heights

  widths.all <- do.call(unit.c, wg) # all grob widths
  heights.all <- do.call(unit.c, hg)    #all grob heights
 
  widths <- colMax.units(widths.all, ncol) + padding.h  # all column widths
  heights <- rowMax.units(heights.all, nrow) + padding.v # all row heights

  if(equal.width)                      
    widths <- rep(max(widths), length(widths))
  if(equal.height)
    heights <- rep(max(heights), length(heights))

  gcells = frameGrob(name="table.cells", vp = "cells",
    layout = grid.layout(nrow, ncol, width=widths, height=heights, just=just) )

  label.ind <- 1   # index running accross labels

  for (ii in seq(1, ncol, 1)) {
    for (jj in seq(1, nrow, 1)) {

      gcells = placeGrob(gcells, lf[[label.ind]], row=jj, col=ii)
      gcells = placeGrob(gcells, lg[[label.ind]], row=jj, col=ii)

      label.ind <- label.ind + 1
    }
  }

  gl = gList( gcells)

  gl
}

exprGrob <- function(e, ncol=n2mfrow(length(e))[2], nrow=n2mfrow(length(e))[1], ...) {
   gTree(e=e, ncol=ncol, nrow=nrow, cl="expr", ...)
}

drawDetails.expr <- function(x, recording=TRUE){
grid.draw(gTree(children=makeTableGrobs(x$e, x$ncol, x$nrow),
                childrenvp=viewport(name="cells")))

}

grid.expr <- function(...){
  grid.draw(exprGrob(...))

}
