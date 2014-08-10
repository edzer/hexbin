\name{panel.hexgrid}
\alias{panel.hexgrid}
\title{Hexagonal grid for a lattice plot}
\description{
A panel function to add a hexagonal grid to a lattice plot.
}
\usage{
panel.hexgrid(h, border = grey(0.85))
}

\arguments{
  \item{h}{an object of class \code{hexbin}.}
  \item{border}{a color for the hexagon border colors}
}
\value{
 There is no return value from this function. The results are plotted on
 the current active device.
}
\author{Nicholas Lewin-Koh \email{nikko@hailmail.net}}
\seealso{\code{\link{hexbinplot}}, \code{\link{hexGraphPaper}}}
\keyword{hplot}
