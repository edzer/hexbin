\name{hexTapply}
\alias{hexTapply}

\title{Apply function to data from each hexagon bin.}

\description{
  A wrapper for tapply except that it operates with each hexagon bin
  being the category. The function operates on the data associated on
  the points from each bin. 
}

\usage{
hexTapply(hbin, dat, FUN = sum, ..., simplify=TRUE)
}
\arguments{
  \item{hbin}{a object of class \code{"hexbin"}, typically produced by
    \code{\link{hexbin}(*)}.}
  \item{dat}{A vector of data the same length as \code{hbin@cID}}
  \item{FUN}{the function to be applied.  In the case of functions like
    \code{+}, \code{\%*\%}, etc., the function name must be quoted.  If
    \code{FUN} is \code{NULL}, tapply returns a vector which can be used
    to subscript the multi-way array \code{tapply} normally produces.}
  \item{\dots}{optional arguments to \code{FUN}.}
  \item{simplify}{If \code{FALSE}, \code{tapply} always returns an array
    of mode \code{"list"}.  If \code{TRUE} (the default), then if
    \code{FUN} always returns a scalar, \code{tapply} returns an array
    with the mode of the scalar.}
}
\details{
  This function is a wrapper for tapply, except that the cell id is
  always the categorical variable. This function is specifically good for
  adding variables to the cAtt slot of a hexbin object or for plotting
  a third variable in a hexagon plot. See below for examples.
}

\value{
  Returns a vector of the result of 'FUN' as in
  \code{\link{tapply}}. See \code{\link{tapply}} for detailed
  description of output.
}

\author{Nicholas Lewin-Koh}
\seealso{ \code{\link{tapply}},\code{\link{hexbin}} }
\examples{
  data(NHANES)
  hbin<-hexbin(log(NHANES$Diet.Iron+1),log(NHANES$BMI),xbins=25,IDs=TRUE)
  hvp<-plot(hbin)
  mtrans<-hexTapply(hbin,NHANES$Transferin,median,na.rm=TRUE)
  pushHexport(hvp$plot.vp)
  grid.hexagons(hbin,style='lattice',pen=0,border='red',use.count=FALSE,
cell.at=mtrans)

  
  
}
\keyword{dplot}
\keyword{utilities}% at least one, from doc/KEYWORDS

