\name{hcell2xyInt}
\alias{hcell2xyInt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Change cell ids to 2d integer coordinate system}
\description{
  Transforms the cell representation of a a lattice into a 2d integer
  coordinate system.
}
\usage{
  hcell2xyInt(hbin, xbins=NULL, xbnds=NULL, ybnds=NULL, shape=NULL)
}

\arguments{
  \item{hbin}{a object of class \code{"hexbin"}, typically produced by
    \code{\link{hexbin}(*)}.}
  \item{xbins}{the number of bins partitioning the range of xbnds.}
  \item{xbnds, ybnds}{horizontal and vertical limits of the binning
    region in x or y units respectively; must be numeric vector of
    length 2.}
  \item{shape}{the \emph{shape} = yheight/xwidth of the plotting regions.}

}
\details{
  Takes a grid defined by either the hexbin parameters or dimen in a
  hexbin object and translates the cell ids for the grid into 2d integer
  coordinates. 
}
\value{
  An integer matrix with two columns, i and j representing the integer
  xy coordinates of the hexagon grid.
  \item{i}{Integer coordiante of the rows, increases from bottom to top}
  \item{j}{Integer coordiante of the columns, increases from left to right}
}

\author{Nicholas Lewin-Koh }


\seealso{\code{\link{hcell2xy}}}
\examples{
 x<-rnorm(10000)
 y<-rnorm(10000)
 hbin<-hexbin(x,y)
 ijInt<-hcell2xyInt(hbin)
}

\keyword{dplot}
\keyword{misc}
