\name{inout.hex}
\alias{inout.hex}
\title{Check points for inclusion}
\description{
  Check which points are in hexagons with \code{count} <= mincnt.
}
\usage{
inout.hex(hbin, mincnt)
}

\arguments{
  \item{hbin}{an object of class \code{\link{hexbin}}.}
  \item{mincnt}{Cutoff, id's for counts less than mincnt are returned}
}
\details{
  Check which points are in hexagons with \code{count} <= mincnt and
  returns the row ids for those points. One can use the ids to plot low
  ount hexagons as points instead.  
}
\value{
  A vector with the row ids of points which fall in hexagons with
  \code{count} less than or equal to mincnt
}

\author{Nicholas Lewin-Koh}


\seealso{\code{\link{plotMAhex}}}

\keyword{misc}

