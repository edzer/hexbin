\name{panel.hexloess}
\alias{panel.hexloess}
\title{Loess line for hexbin lattice plot}
\description{
  A panel function to add a loess line to a hexbin lattice plot.

  This function contravened CRAN policy and is no longer available.
}
\usage{
panel.hexloess(bin, w = NULL, span = 2/3, degree = 1,
               family = c("symmetric", "gaussian"), evaluation = 50, 
               lwd = add.line$lwd, lty = add.line$lty,
               col, col.line = add.line$col, \dots)

}
\arguments{
  \item{bin}{an object of class \code{hexbin}.}
  \item{w}{optional counts for object \code{bin}.}
  \item{span}{smoothness parameter for \code{loess}.}
  \item{degree}{degree of local polynomial used.}
  \item{family}{if \code{"gaussian"} fitting is by least-squares, and
                if \code{"symmetric"} a re-descending M-estimator is used.}
  \item{evaluation}{number of points at which to evaluate the smooth curve.}
  \item{lwd}{line weight graphical parameter.}
  \item{lty}{line type graphical parameter.}
  \item{col}{same as \code{col.line}.}
  \item{col.line}{line color graphical parameter.}
  \item{\dots}{optional arguments to \code{\link[stats]{loess.control}}.}
}
\value{
 There is no return value from this function. The results are plotted on
 the current active device.
}
\author{Nicholas Lewin-Koh \email{nikko@hailmail.net}}
\seealso{
\code{\link{hexbinplot}}, \code{\link{panel.hexgrid}},
\code{\link{loess.smooth}},
\code{\link{loess.control}}, 
\code{\link[lattice:panel.functions]{panel.loess}}
}
\keyword{hplot}
