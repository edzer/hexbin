\name{hexpolygon}
\alias{hexpolygon}
\alias{hexcoords}
\title{Hexagon Coordinates and Polygon Drawing}
\description{
  Simple \sQuote{low-level} function for computing and drawing hexagons.
  Can be used for \sQuote{grid} (package \pkg{grid}) or
  \sQuote{traditional} (package \pkg{graphics}) graphics.
}
\usage{
hexcoords(dx, dy = NULL, n = 1, sep = NULL)

hexpolygon(x, y, hexC = hexcoords(dx, dy, n = 1), dx, dy = NULL,
           fill = 1, border = 0, hUnit = "native", ...)
}
\arguments{
  \item{dx,dy}{horizontal and vertical width of the hexagon(s).}
  \item{n}{number of hexagon \dQuote{repeats}.}
  \item{sep}{separator value to be put between coordinates of different
    hexagons. The default, \code{NULL} doesn't use a separator.}
  \item{x,y}{numeric vectors of the same length specifying the hexagon
    \emph{centers} around which to draw.}
  \item{hexC}{a list as returned from \code{hexcoords()}.
    Its component \code{no.sep} determines if grid or traditional
    graphics are used. The default (via default of \code{hexcoords}) is
    now to use grid graphics.}
  \item{fill,border}{passed to \code{\link[grid]{grid.polygon}} (for \pkg{grid}).}
  \item{hUnit}{string or \code{NULL} determining in which units (x,y)
    values are.}
  \item{\dots}{further arguments passed to \code{\link{polygon}} (for
    \pkg{graphics}).}
}
\value{
  \code{hexcoords()} returns a list with components
  \item{x,y}{numeric vectors of length \eqn{n \times 6}{n * 6} (or
    \eqn{n \times 7}{n * 7} if \code{sep} is not NULL)
    specifying the hexagon polygon coordinates (with \code{sep} appended
    to each 6-tuple).}
  \item{no.sep}{a logical indicating if \code{sep} was \code{NULL}.}

  \code{hexpolygon} returns what its last \code{\link[grid]{grid.polygon}(.)}
  or \code{\link{polygon}(.)} call returns.
}
\author{Martin Maechler, originally.}
\seealso{\code{\link{grid.hexagons}} which builds on these.}
\examples{
str(hexcoords(1, sep = NA))  # multiple of (6 + 1)
str(hexcoords(1, sep = NULL))# no separator -> multiple of 6
\dontshow{
stopifnot(3 * (6+1) == sapply(hexcoords(2, n = 3, sep = NA)[1:2], length),
          6 == sapply(hexcoords(1)[1:2], length))
set.seed(1001)
}

## hexpolygon()s:
x <- runif(20, -2, 2)
y <- x + rnorm(20)

## 1) traditional 'graphics'
plot(x,y, asp = 1, "plot() + hexpolygon()")
hexpolygon(x,y, dx = 0.1, density = 25, col = 2, lwd = 1.5)

## 2) "grid" :

addBit <- function(bnds, f = 0.05) bnds + c(-f, f) * diff(bnds)
sc <- addBit(rxy <- range(x,y))# same extents (cheating asp=1)
library("grid")
grid.newpage()
pushViewport(plotViewport(.1+c(4,4,2,1), xscale = sc, yscale = sc))
grid.rect()
grid.xaxis()
grid.yaxis()
grid.points(x,y)
hexpolygon(x,y, hexcoords(dx = 0.1, sep=NULL), border = "blue", fill=NA)
popViewport()
}
\keyword{dplot}
\keyword{aplot}
