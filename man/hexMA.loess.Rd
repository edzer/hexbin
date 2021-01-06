\name{hexMA.loess}
\alias{hexVP.loess}
\alias{hexMA.loess}
\title{Add Loess Fit to Hexplot }
\description{
  Fit a loess line using the hexagon centers of mass as the x and y
  coordinates and the cell counts as weights.
}
\usage{
hexMA.loess(pMA, span = 0.4, col = "red", n = 200, ...)
hexVP.loess(hbin, hvp = NULL, span = 0.4, col = "red", n = 200, ...)
}

\arguments{
  \item{hbin}{an object of class \code{hexbin}, see \code{\link{hexbin}}.}
  \item{hvp}{A \code{hexViewport} object.}
  \item{pMA}{the list returned by \code{\link{plotMAhex}}.}
  \item{span}{the parameter alpha which controls the degree of smoothing.}
  \item{col}{line color for the loess fit.}
  \item{n}{number of points at which the fit should be evaluated.}
  \item{...}{Additional graphical parameter settings for the \code{loess} line fit; see \code{\link[grid]{gpar}}.}
}
  
\value{
  Returns invisibly the object associated with the loess fit.
}

\author{Nicholas Lewin-Koh }

\seealso{ \code{\link{hexVP.abline}}, \code{\link{plotMAhex}},
  \code{\link{gplot.hexbin}}, \code{\link{hexViewport}};
  \code{\link{loess}}
}
\examples{
  if(require(marray)){
    data(swirl)
    %% the following had     'newpage = FALSE, '  -- why ??
    hb <- plotMAhex(swirl[,1], main = "M vs A plot with hexagons", legend=0)
    hexVP.abline(hb$plot, h=0, col= gray(.6))
    hexMA.loess(hb)
  }
  dat <- data.frame(x=rnorm(1000), y=rnorm(1000))
  bin <- hexbin(dat$x, dat$y)
  hb <- plot(bin)
  hexVP.loess(bin, hvp = hb$plot.vp, span = 0.4, n = 200, col = "blue", lwd = 3, lty = "dashed")
}
\keyword{aplot}

