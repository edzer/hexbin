\name{hexVP.abline}
\alias{hexVP.abline}
\title{Add a Straight Line to a HexPlot}
\description{
  This function adds one or more straight lines through the current
  plot; it is the hexbin version of \code{\link[graphics]{abline}()}.
}
\usage{
hexVP.abline(hvp, a = NULL, b = NULL, h = numeric(0), v = numeric(0),
             col = "black", lty = 1, lwd = 2, \dots)
}
\arguments{
  \item{hvp}{A hexViewport object that is currently on the active device}
  \item{a,b}{the intercept and slope or if \code{b} is \code{NULL},
    an \code{lm} object or a vector of length 2 with
    \code{c(intercept,slope)}}
  \item{h}{the y-value for a horizontal line.}
  \item{v}{the x-value for a vertical line.}
  \item{col, lty, lwd}{line color, type and width.}
  \item{\dots}{further graphical parameters.}
}
\details{
  The first form specifies the line in intercept/slope form
  (alternatively \code{a} can be specified on its own and is taken to
  contain the slope and intercept in vector form).

  The \code{h=} and \code{v=} forms draw horizontal and vertical
  lines at the specified coordinates.

  The \code{coef} form specifies the line by a vector containing the
  slope and intercept.

  \code{lm} is a regression object which contains \code{reg$coef}.  If it is
  of length 1 then the value is taken to be the slope of a line
  through the origin, otherwise, the first 2 values are taken to be
  the intercept and slope.
}
\author{Nicholas Lewin-Koh}

\seealso{\code{\link{gplot.hexbin}}, \code{\link{hexViewport}},
  \code{\link{hexMA.loess}}
}
\keyword{aplot}

