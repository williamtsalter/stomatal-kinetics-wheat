library(grid)

ceiling2 <- function(x,s=1)  s*ceiling(x/s)
floor2   <- function(x,s=1)  s*floor(x/s)
mod.pos  <- function(x,s) (x-1)%%s + 1

rcpos <- function(row=1,col=1)  viewport(layout.pos.row=row, layout.pos.col=col)
rpos  <- function(row=1)  rcpos(row,1)
cpos  <- function(col=1)  rcpos(1,col)

makeGrid <- function(rows=1,cols=1) {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(rows,cols)))
}

plotPanel <- function(plt, rows, cols, panel=1, dir="v")
{
  if (dir=="v") # direction = vertical
    print(plt, vp = rcpos(row = mod.pos(panel,rows), col = ceiling(panel/rows)))
  else          # direction = horizontal
    print(plt, vp = rcpos(col = mod.pos(panel,cols), row = ceiling(panel/cols)))
}

# Useful in Rmd files when multiple figures are created in a single code-chunk
drawFooter <- function()  cat("##  ---  ##")


doubleSummed <- function(tbl) {
  tbl <- cbind(tbl, rowSums(tbl))
  tbl <- rbind(tbl, colSums(tbl))
  colnames(tbl)[ncol(tbl)] <- "ALL"
  rownames(tbl)[nrow(tbl)] <- "ALL"
  tbl
}

doubleSummary <- function(data, y, x1, x2, FUN=mean, ...)
{
  yval <- data[[y]]
  rval <- data[[x1]]
  cval <- data[[x2]]
  
  tbl <- rbind(
    cbind(
      tapply(yval, list(rval,cval), FUN, ...),
      tapply(yval, list(rval), FUN, ...)),
    c(tapply(yval, list(cval), FUN, ...),
      FUN(yval, ...))
  )
  colnames(tbl)[ncol(tbl)] <- "ALL"
  rownames(tbl)[nrow(tbl)] <- "ALL"
  tbl
}

addRates <- function(resp,sep=".")  c(resp, paste(resp,c("AGR","RGR"),sep=sep))

intervalTraits.str <- function(intervals, responses, sep)
{
  as.vector(outer(intervals, responses, function(x,y) paste(y,x,sep=sep)))
}

intervalTraits <- function(starts, ends, response="Area.smooth", sep1="to", sep2=".")
{
  intervalTraits.str(paste(starts, ends, sep=sep1),
                     paste(response, c("AGR","RGR"), sep=sep2),
                     sep = sep2)
}


make.scale <- function(zmin,zmax,dz1,dz2)
{
  dz2 = min(dz2,dz1)
  
  val <- list()
  val$limits <- c(zmin,zmax)
  val$labels <- seq(ceiling2(zmin,dz1), zmax, by=dz1)
  val$lines  <- seq(ceiling2(zmin,dz2), zmax, by=dz2)
  
  val
}

make.xscale <- function(xmin, xmax, dx1, dx2=dx1, maxBreaks=51, limit.breaks=TRUE, ...)
{
  val <- make.scale(xmin,xmax,dx1,dx2)
  
  if (limit.breaks & length(val$labels) > maxBreaks) {
    warning("Horizontal grid spacing dx1 too small, reverting to automatic scale.")
    return(NULL)
  }
  
  scale_x_continuous(limits = val$limits,
                     breaks = val$labels,
                     minor_breaks = val$lines, ...)
}

make.yscale <- function(ymin, ymax, dy1, dy2=dy1, maxBreaks=51, limit.breaks=TRUE, ...)
{
  val <- make.scale(ymin,ymax,dy1,dy2)
  
  if (limit.breaks & length(val$labels) > maxBreaks) {
    warning("Vertical grid spacing dy1 too small, reverting to automatic scale.")
    return(NULL)
  }
  
  scale_y_continuous(limits = val$limits,
                     breaks = val$labels,
                     minor_breaks = val$lines, ...)
}

data.xscale <- function(data, dx1, dx2=dx1,
                        xmin =   floor2(min(data, na.rm=T), dx2),
                        xmax = ceiling2(max(data, na.rm=T), dx2), 
                        maxBreaks=51, limit.breaks=TRUE, ...)
{
  make.xscale(xmin, xmax, dx1, dx2, maxBreaks=maxBreaks, limit.breaks=limit.breaks, ...)
}

data.yscale <- function(data, dy1, dy2=dy1,
                        ymin =   floor2(min(data, na.rm=T), dy2),
                        ymax = ceiling2(max(data, na.rm=T), dy2), 
                        maxBreaks=51, limit.breaks=TRUE, ...)
{
  make.yscale(ymin, ymax, dy1, dy2, maxBreaks=maxBreaks, limit.breaks=limit.breaks, ...)
}

make.ylogscale <- function(ymin, ymax, constrained=FALSE, coarse=FALSE, ...)
{
  a  <-   floor(log10(ymin))
  b  <- ceiling(log10(ymax))
  xx <- 10^(a:(b-1))
  
  if (coarse) {
    major <- 10^(a:b)
    minor <- as.numeric(2:5 %o% xx)
  } else {
    major <- as.numeric(c(c(1,2,5)  %o% xx, 10^b))
    minor <- as.numeric(c(1.5, 3:9) %o% xx)
  }
  
  if (constrained)
    scale_y_log10(breaks=major, minor_breaks=minor, limits=c(ymin,ymax), ...)
  else
    scale_y_log10(breaks=major, minor_breaks=minor, ...)
}

make.xlogscale <- function(xmin, xmax, constrained=FALSE, coarse=FALSE, ...)
{
  a  <-   floor(log10(xmin))
  b  <- ceiling(log10(xmax))
  xx <- 10^(a:(b-1))
  
  if (coarse) {
    major <- 10^(a:b)
    minor <- as.numeric(2:5 %o% xx)
  } else {
    major <- as.numeric(c(c(1,2,5)  %o% xx, 10^b))
    minor <- as.numeric(c(1.5, 3:9) %o% xx)
  }
  
  if (constrained)
    scale_x_log10(breaks=major, minor_breaks=minor, limits=c(xmin,xmax), ...)
  else
    scale_x_log10(breaks=major, minor_breaks=minor, ...)
}
