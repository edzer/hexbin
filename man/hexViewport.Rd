\name{hexViewport}
\alias{hexViewport}
\title{Compute a Grid Viewport for Hexagon / Hexbin Graphics}
\description{
  Builds a \code{grid} viewport for hexagon or \code{\link{hexbin}}
  graphics.  This builds on the concepts of the \pkg{grid} package,
  see \code{\link[grid]{viewport}}.
}
\usage{% see ../R/hexViewport.R
hexViewport(x, offset = unit(0,"inches"), mar = NULL,
	    xbnds = NULL, ybnds = NULL, newpage = FALSE,
            clip = "off", vp.name = NULL)
}
\arguments{
  \item{x}{a \code{\link{hexbin}} object.}
  \item{offset}{a \code{\link[grid]{unit}} object.}
  \item{mar}{margins as \code{\link[grid]{unit}}s, of length 4 or 1.}
  \item{xbnds, ybnds}{bounds for x- and y- plotting range; these default
    to the corresponding slots of \code{x}.}
  \item{newpage}{logical indicating if a new graphics page should be
    openend, i.e., \code{\link[grid]{grid.newpage}()}.}
  \item{clip}{simply passed to \code{\link[grid]{viewport}()}.}
  \item{vp.name}{name of viewport; defaults to random name.}
}
\value{
  an S4 object of class \code{"hexVP"}, see \link{hexVP-class} for more,
  with its main slot \code{hexVp} a \code{\link[grid]{viewport}} for
  grid graphics. 
}

\seealso{\code{\link[grid]{viewport}} and the main
  \emph{\dQuote{handlers}} \code{\link{pushHexport}} and
  \code{\link[grid:viewports]{popViewport}}; further
  \code{\link{gplot.hexbin}} and \code{\link{hboxplot}} which build on
  \code{hexViewport}.
}
\examples{
set.seed(131)
x <- rnorm(7777)
y <- rt   (7777, df=3)

## lower resolution binning and overplotting with counts
bin <- hexbin(x,y,xbins=25)
P <- plot(bin)
xy <- hcell2xy(bin)
pushHexport(P$plot.vp)
i <- bin@count <= 3
library("grid")
grid.text(as.character(bin@count[i]), xy$x[i], xy$y[i],
          default.units = "native")
grid.points(x[1:20],y[1:20]) # to show some points rather than counts
popViewport()
}
\keyword{hplot}% ?
\keyword{aplot}
