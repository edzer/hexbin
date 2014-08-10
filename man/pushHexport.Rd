\name{pushHexport}
\alias{pushHexport}
\title{Push a Hexagon Viewport ("hexVP")}
\description{
  Push a Hexagon Viewport (\code{"hexVP"}, see \link{hexVP-class}) on to
  the tree of (grid) viewports, calling
  \code{\link[grid:viewports]{pushViewport}}. 
}
\usage{
pushHexport(hvp, clip = "off")
}
\arguments{
  \item{hvp}{a hexagon viewport, i.e., an object of class
    \code{"hexVP"}, see \link{hexVP-class}, typically produced by
    \code{\link{hexViewport}(..)}.}
  \item{clip}{which viewport to push, either 'on' or 'off' are the
    allowed arguments, see details.}
}
\seealso{the underlying \code{\link[grid:viewports]{pushViewport}} from the
  \pkg{grid} package.
}
\details{
  A hexagon viewport (\code{"hexVP"}) object has slots for two replicate
  viewports one with clipping turned on and one with clipping off.  This
  allows toggling the clipping option.
}
%\examples{
%}
\keyword{dplot}
