\name{getHMedian}
\alias{getHMedian}
\alias{getHMedian,erodebin-method}
\title{Get coordiantes of the median cell after the erode operation}
\description{
  A method for a eroded hexbin object to extract the coordinates of the
  median cell.  The median is simply the cell with the highest erosion
  number or the last cell to be eroded.
}
\usage{
getHMedian(ebin)
}
\arguments{
  \item{ebin}{result of \code{\link{erode.hexbin}()}.}
}
\section{Methods}{
  \describe{
    \item{ebin = "erodebin"}{...}
  }
}
\seealso{\code{\link{erode.hexbin}}
}
\examples{
set.seed(153)
x <- rnorm(10000)
y <- rnorm(10000)
bin <- hexbin(x,y)

smbin  <- smooth.hexbin(bin)
erodebin <- erode.hexbin(smbin, cdfcut=.5)
getHMedian(erodebin)
}
\keyword{methods}

