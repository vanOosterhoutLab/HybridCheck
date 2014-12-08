#' display a data.frame
#'
#' uses Grid graphics to display a data.frame on a graphics page
#' @aliases tableGrob2 grid.table2
#' @title tableGrob2
#' @param d data.frame
#' @param name  
#' @param gp  
#' @param show.rownames logical
#' @param show.colnames logical
#' @param equal.width logical
#' @param core.vsep logical 
#' @param core.hsep logical
#' @param head.vsep logical 
#' @param head.hsep logical 	
#' @param core.fill colour
#' @param head.fill colour
#' @param core.border colour
#' @param head.border colour
#' @param core.text colour
#' @param head.text colour
#' @return a grob
#' @seealso \code{grid.segments}, \code{grid.points} 
#' @examples
#' tc  = textConnection("
#'      carat   VeryLongWordIndeed color clarity depth
#' 14513  1.35 Ideal     J     VS2  61.4
#' 28685  0.30  Good     G    VVS1  64.0
#' 50368  0.75 Ideal     F     SI2  59.2")
#' d = read.table(tc,head=T)
#' close(tc)
#' grid.newpage()
#' g = grid.table2(d)
#' grid.ls(g) 
#' grid.edit("top-head-fill-5", gp=gpar(fill="red"))
#' grid.edit("cells-label-33", label=expression(alpha),gp=gpar(col="orange"))
#' grid.newpage()
#' grid.table2(d, core.fill ="red", head.text = "green",
#' 	core.text = "blue", head.fill ="cyan",head.border="black",core.border="magenta")
#' grid.newpage()
#' grid.table2(d,equal=T)
#' grid.newpage()
#' grid.table2(d,show.rownames=F)
#' grid.newpage()
#' grid.table2(d,show.colnames=F)


makeTableGrobs2 <- function(d, widths, left.width, content.width, top.height, content.height,
	just = c("left", "top"),
	show.rownames=T,show.colnames=T,
		core.vsep=T, core.hsep=T, 	
		head.vsep=T, head.hsep=T, 	
		core.fill = "grey95",
			head.fill = "grey90",
			core.border = "white",
			head.border = "white",
			core.text = "black",
			head.text = "black") {
			
	ncol = ncol(d)
	nrow = nrow(d)
	
# inner cells	
	gcells = frameGrob(name="table.cells", vp = "cells", 
		layout = grid.layout(nrow, ncol, just=just,
			widths = widths,
			heights = unit(rep(1, nrow), "lines") ))
			
	for (ii in seq(1, ncol, 1)) {
		for (jj in seq(1, nrow, 1)) {
		gcells = placeGrob(gcells, rectGrob(gp=gpar(fill = core.fill, col=NA), name=paste("cells-fill-",ii,jj,sep="")), row=jj, col=ii)
		if (core.hsep )
		gcells = placeGrob(gcells, segmentsGrob(0,0,0,1, gp=gpar(col = core.border)), 
			row=jj, col=ii)
		if (core.vsep)
		gcells = placeGrob(gcells, segmentsGrob(0,0,1,0, gp=gpar(col = core.border)), 
			row=jj, col=ii)
		gcells = placeGrob(gcells, textGrob(label=d[jj,ii], gp=gpar(col=core.text), name=paste("cells-label-",ii,jj,sep="")), 
			row=jj, col=ii)
		}
	}

if (show.colnames)
{
# top header
	gtop = frameGrob(name="table.header.top", vp = "header.top",
	   	layout = grid.layout(1, ncol, just=just,
	   		widths = widths,
	   		heights = unit(1, "lines") ))
	
 for (ii in seq(1, ncol, 1)) {
 	gtop = placeGrob(gtop, rectGrob(gp=gpar(fill = head.fill,col=NA), name=paste("top-head-fill-",ii,sep="")), row=1, col=ii)
 	if (head.hsep) 
 	gtop = placeGrob(gtop, segmentsGrob(0,0,0,1,gp=gpar(col = head.border)), row=1, col=ii)
 	if (head.vsep) 
 	gtop = placeGrob(gtop, segmentsGrob(0,0,1,0,gp=gpar(col = head.border)), row=1, col=ii)
 	gtop = placeGrob(gtop, 
		textGrob(label=colnames(d)[ii], gp=gpar(col=head.text,fontface="bold"), name=paste("top-head-label-", ii,sep="")), row=1, col=ii)
 }
}


if (show.rownames)
{
# left header
	 gleft = frameGrob(name="table.header.left", vp="header.left",
	 layout = grid.layout(nrow, 1, just=just, 
	 	widths = left.width,
	 	heights = unit(1, "lines") ))

 for (jj in seq(1, nrow, 1)) {
 	gleft = placeGrob(gleft, rectGrob(gp=gpar(fill = head.fill,col=NA), name=paste("left-head-fill-", jj,sep="")), row=jj, col=1)
 	if (head.hsep)
 	gleft = placeGrob(gleft, segmentsGrob(1,0,1,1, gp=gpar(col = head.border)), row=jj, col=1)
 	if (head.vsep)
 	gleft = placeGrob(gleft, segmentsGrob(0,1,1,1, gp=gpar(col = head.border)), row=jj, col=1)
 	gleft = placeGrob(gleft, 
		textGrob(label=rownames(d)[jj],gp=gpar(col=head.text,fontface="bold"), name=paste("left-head-label-", jj,sep="")), row=jj, col=1)
 }
}
				
if (show.rownames & show.colnames)
gl = gList( gcells, gleft, gtop)

if (!show.rownames & !show.colnames)
gl = gList( gcells)

if (show.rownames & !show.colnames)
gl = gList( gcells, gleft)

if (!show.rownames & show.colnames)
gl = gList( gcells, gtop)			

gl
}




makeTableVps2 <- function(show.rownames=T,show.colnames=T) {

if (show.rownames & show.colnames)
vl = vpList(viewport(name="cells", layout.pos.row = 2, layout.pos.col = 2), 
	   viewport(name="header.top", layout.pos.row = 1, layout.pos.col = 2),
	   viewport(name="header.left", layout.pos.row = 2, layout.pos.col = 1)
)

if (!show.rownames & !show.colnames)
vl = vpList(viewport(name="cells", layout.pos.row = 1, layout.pos.col = 1)
)

if (show.rownames & !show.colnames)
vl = vpList(viewport(name="cells", layout.pos.row = 1, layout.pos.col = 2), 
	   viewport(name="header.left", layout.pos.row = 1, layout.pos.col = 1)
)

if (!show.rownames & show.colnames)
vl = vpList(viewport(name="cells", layout.pos.row = 2, layout.pos.col = 1), 
	   viewport(name="header.top", layout.pos.row = 1, layout.pos.col = 1)
)

vl
}
				

tableGrob2 <- function(d,  
                    name=NULL, gp=NULL, show.rownames=T, show.colnames=T, equal.width=F,
					just = c("left", "top"),
					core.vsep=T, core.hsep=T, 	
					head.vsep=T, head.hsep=T, 	
					core.fill = "grey95",
					head.fill = "grey90",
					core.border = "white",
					head.border = "white",
					core.text = "black",
					head.text = "black") { 
						
	ncol = ncol(d)
	nrow = nrow(d)
	
	dm = as.matrix(d) 
	dm = cbind(rownames(d), dm)
	dm = rbind(colnames(dm), dm)

    left.width = convertUnit(max(stringWidth(rownames(dm))), "mm") + rep(unit(4,"mm"), ncol )
    top.height = unit(1,"line")
    content.height = unit(nrow,"line")
	
	if(equal.width){
	widths = rep(convertUnit(max(stringWidth(dm)), "mm") + unit(4,"mm"), ncol )
	left.width = widths[1]
	}
	
	if(!equal.width)
	widths <- unit(4,"mm") + 
	do.call(unit.c, lapply(seq(2, ncol(dm), 1), 
		function(.col) convertUnit(max(stringWidth(dm[seq(1+!show.colnames, nrow(dm), 1),.col])), "mm")))

    content.width = sum(convertUnit(widths, "mm"))
		
	heights = rep( unit(1,"line"), nrow + show.colnames)
	
	ws = if(show.rownames) unit.c(left.width,content.width) else content.width
	hs = if(show.colnames) unit.c(top.height, content.height) else content.height
	vp.table =  viewport(name="table", layout=
		grid.layout(1+show.colnames,1+show.rownames, width = ws, 
						height = hs, just=c("left", "top")
				   )
				)
  
		
  g <- gTree(nrow=nrow, ncol=ncol, just = just, show.rownames=show.rownames, show.colnames=show.colnames, 
               children = makeTableGrobs2(d, widths, left.width, content.width, top.height, content.height,
										just = just,
			  							show.rownames=show.rownames, show.colnames=show.colnames, 
			  							core.vsep=core.vsep, core.hsep=core.hsep, 	
			  							head.vsep=head.vsep, head.hsep=head.hsep, 	
			  							core.fill = core.fill,
			  							head.fill = head.fill,
			  							core.border = core.border,
			  							head.border = head.border,
			  							core.text = core.text,
			  							head.text = head.text),
               gp=gp, name=name, vp = vp.table, 
			   childrenvp = makeTableVps2(show.rownames=show.rownames,show.colnames=show.colnames),
               cl = "tableGrob") 
  g
}

grid.table2 <- function(...) {
  g <- tableGrob2(...)
  grid.draw(g)
invisible(g)
}





