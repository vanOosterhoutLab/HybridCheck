

grobii <- function(d, gp=gpar(), name="row-label-",
                   just="center", parse=TRUE){
  x <- switch(just, "center"=0.5, "right"=0.95, "left"=0.1)
   function(ii){
     g <- d[[ii]]
     gTree(x=x, children=gList(g), just=just, gp=gp, name=paste(name, ii, sep=""))
  }
}

makeTableGrobs <- function(content, rnames=NULL, cnames=NULL,
        nrow, ncol, parse=TRUE,
        row.just="center", col.just="center", core.just="center",
        equal.width = FALSE, equal.height=FALSE, 
        gpar.coretext = gpar(col="black", cex=1),
        gpar.coltext = gpar(col="black", cex=1, fontface="bold"),
        gpar.rowtext = gpar(col="black", cex=0.8, fontface="italic"),
        h.odd.alpha = 1, h.even.alpha = 1, 
        v.odd.alpha = 1, v.even.alpha = 1, 
        gpar.corefill = gpar(fill = "grey95", col="white"), 
        gpar.rowfill = gpar(fill = "grey90", col="white"), 
        gpar.colfill = gpar(fill = "grey90", col="white")) {

 ncontent <- length(content) # number of labels
 nrnames <- length(rnames) # number of row labels
 ncnames <- length(cnames) # number of col labels

 
 ## define some functions to generate named grobs
 makeOneRowname <- if(all(is.character(rnames)))
   textii(d=rnames, gp=gpar.rowtext, name="row-label-", just=row.just, parse=parse) else 
 grobii(d=rnames, gp=gpar.rowtext, name="row-label-", just=row.just)
 
 makeOneColname <- if(all(is.character(cnames)))
   textii(d=cnames, gp=gpar.coltext, name="col-label-", just=col.just, parse=parse) else
 grobii(d=cnames, gp=gpar.coltext, name="col-label-", just=col.just)
 
 makeOneLabel <- textii(d=content, gp=gpar.coretext, name="core-label-", just=core.just, parse=parse)
 

 gp.corefillee <- gp.corefilleo <- gp.corefilloe <- gp.corefilloo <- gpar.corefill
 gp.corefillee[["alpha"]] <- h.even.alpha *  v.even.alpha
 gp.corefilloe[["alpha"]] <- h.odd.alpha *  v.even.alpha
 gp.corefilloo[["alpha"]] <- h.odd.alpha *  v.odd.alpha
 gp.corefilleo[["alpha"]] <- h.even.alpha *  v.odd.alpha
 
 gpar.corefill <- rep(c(rep(c(list(gp.corefillee), list(gp.corefilloe)), length.out=nrow),
                        rep(c(list(gp.corefilleo), list(gp.corefilloo)), length.out=nrow)),
                      length.out=ncontent)
 
 gp.rowfille <- gp.rowfillo <- gpar.rowfill
 gp.rowfille[["alpha"]] <-  h.even.alpha
 gp.rowfillo[["alpha"]] <-  h.odd.alpha

 gpar.rowfill <- rep(c(list(gp.rowfille), list(gp.rowfillo)), nrow)
 
 gp.colfille <- gp.colfillo <- gpar.colfill
 gp.colfille[["alpha"]] <-  v.even.alpha
 gp.colfillo[["alpha"]] <-  v.odd.alpha
 
 gpar.colfill <- rep(c(list(gp.colfille), list(gp.colfillo)), ncol)
   
 makeOneCell <- rectii(gp=gpar.corefill, name="core-fill-")
 makeOneRowfill <- rectii(gp=gpar.rowfill, name="row-fill-")
 makeOneColfill <- rectii(gp=gpar.colfill, name="col-fill-")

 ## in case of missing row(col) names,  make a list of virtualGrobs
 ## else,  a list of rectGrobs with incremental names
 if(is.null(rnames)){
   lrt <- lrf <- replicate(nrow, virtualGrob, simplify=FALSE)} else {
  lrt <- lapply(seq_along(rnames), makeOneRowname) # list of text grobs
  lrf <- lapply(seq_along(rnames), makeOneRowfill) # list of rect grobs
   }
 if(is.null(cnames)){
  lct <- lcf <- replicate(ncol, virtualGrob, simplify=FALSE)} else {
  lct <- lapply(seq_along(cnames), makeOneColname) # list of text grobs
  lcf <- lapply(seq_along(cnames), makeOneColfill) # list of rect grobs
   }
 ## the content consists of textGrobs and rectGrobs
  lit <- lapply(seq_along(content), makeOneLabel) # list of text grobs
  lif <- lapply(seq_along(content), makeOneCell) # list of rect grobs

 ## here the grobs are arranged and permuted in a list to fill a matrix column by column
 lgt <- c(list(virtualGrob), lrt, interleaven(lct, lit, nrow)) # all labels in order
 lgf <- c(list(virtualGrob), lrf, interleaven(lcf, lif, nrow)) # all labels in order
 
 ## retrieve the widths and heights of all textGrobs (including some virtualGrobs) 
  wg <- lapply(lgt, grobWidth) # list of grob widths
  hg <- lapply(lgt, grobHeight) # list of grob heights

 ## concatenate this units
  widths.all <- do.call(unit.c, wg) # all grob widths
  heights.all <- do.call(unit.c, hg)    #all grob heights

 ## matrix-like operations on units to define the table layout
  widths <- colMax.units(widths.all, ncol+1)  # all column widths
  heights <- rowMax.units(heights.all, nrow+1)  # all row heights

 ## equal width or equal height (all except rows and cols)
 nwidths <- length(widths)
 nheights <- length(heights)
  if(equal.width)                      
    widths <- unit.c(widths[[1]], rep(max(widths[seq(2, nwidths)]), nwidths-1))
  if(equal.height)
    heights <- unit.c(heights[[1]], rep(max(heights[seq(2, nheights)]), nheights-1))

 ## return a list containing lists of grobs, and the dimensions for a rectangular layout
 list(lgt=lgt, lgf=lgf, nrow=nrow, ncol=ncol, widths=widths, heights=heights)
}

utils::assignInNamespace("makeTableGrobs", makeTableGrobs, "gridExtra")
