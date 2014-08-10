\name{optShape}
\alias{optShape}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Optimal Shape Parameter for Hexbin Viewport}
\description{
  Takes a viewport or a given height and width and returns the shape
  parameter that will fill the specified plotting region with the
  appropriately shaped hexagons.  If margins are specified the margins
  are subtracted from height and width before the shape parameter is
  specified.
}
\usage{
optShape(vp, height = NULL, width = NULL, mar = NULL)
}
\arguments{
  \item{vp}{a \code{viewport} object, optional see details}
  \item{height}{the height of the plotting region, can be numeric or units}
  \item{width}{The width of the plotting region, can be numeric or units}
  \item{mar}{A four element numeric or units vector describing the
    margins in the order \code{c(bottom, left, top, right)}}
}
\value{
  a scalar numeric value specifiyng \code{shape}.
}
\author{Nicholas Lewin-Koh}
\section{Warning}{If a viewport is given as an argument it should
  already be pushed on the graphics device or it will have null units
  and a meaningless shape parameter will be returned.
}
\seealso{\code{\link{hexViewport}}, \code{\link{hexVP-class}},
         \code{\link{hexbin}}}
\examples{
x <- rgamma(10000,.9)
m <- as.logical(rbinom(10000,1,.17))
x[m] <- -x[m]
y <- rnorm(x,abs(x))
library("grid")
vp <- plotViewport(xscale= range(x)+c(-.5,.5),
                   yscale= range(y)+c(-.5,.5),
             default.units = "native")
grid.newpage()
pushViewport(vp)
grid.rect()
shape <- optShape(vp)
shape
hb <- hexbin(x,y,xbins=40,shape=shape)
grid.hexagons(hb,colramp=BTY)
}
\keyword{dplot}

