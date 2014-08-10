\name{hexList}
\alias{hexList}
\alias{hexbinList-class}
\alias{coerce,list,hexbinList-method}
\title{Conditional Bivariate Binning into Hexagon Cells }
\description{
  Creates a list of \code{\link{hexbin}} objects.  Basic components are
  a cell id and a count of points falling in each occupied cell.
  Basic methods are \code{\link[methods]{show}()}, \code{plot()} %(\link{plot.hexbin})
  and \code{\link{summary}()}, but also \code{\link{erode}}.
  % .. \code{\link{smooth.hexbin}}
}
\usage{
hexList(x, y = NULL, given = NULL, xbins = 30, shape = 1,
        xbnds = NULL, ybnds = NULL, xlab = NULL, ylab = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ x coordinate to be binned }
  \item{y}{ y coordinate to be binned }
  \item{given}{ ..}
  \item{xbins}{ number of bins partitioning the range of xbnds}
  \item{shape}{ the \emph{shape} = yheight/xwidth of the plotting regions }
  \item{xbnds}{ horizontal limits of binning }
  \item{ybnds}{ vertical limits of binning }
  \item{xlab}{ character strings used as labels for \code{x} }
  \item{ylab}{ character strings used as labels for \code{y}}
}
\details{
  There is also a \code{\link[methods:as]{coerce}} method to produce
  \code{hexbinList} objects from \code{\link{list}}s.
  %% i.e., \code{as(list(....), "hexbinList")} should work
}
\value{

  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}

\author{Nicholas Lewin-Koh}

\seealso{\code{\link{hexbin}}, \code{\link{hdiffplot}} }
\keyword{dplot}
\keyword{misc}
