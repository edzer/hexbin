\name{hexbin}
\title{Bivariate Binning into Hexagon Cells}
\alias{hexbin}
\alias{hexbin-class}
\alias{integer or NULL-class}
\alias{show,hexbin-method}
\alias{summary,hexbin-method}
\description{
  Creates a \code{"hexbin"} object.  Basic components are a cell id and
  a count of points falling in each occupied cell.

  Basic methods are \code{\link[methods]{show}()}, \code{plot()} %(\link{plot.hexbin})
  and \code{\link{summary}()}, but also \code{\link{erode}}.
  % .. \code{\link{smooth.hexbin}}
}
\usage{
hexbin(x, y, xbins = 30, shape = 1,
       xbnds = range(x), ybnds = range(y),
       xlab = NULL, ylab = NULL, IDs = FALSE)
}
\arguments{
  \item{x, y}{vectors giving the coordinates of the bivariate data
    points to be binned.  Alternatively a single plotting structure can
    be specified: see \code{\link[grDevices]{xy.coords}}.  \code{\link{NA}}'s are
    allowed and silently omitted.}
  \item{xbins}{the number of bins partitioning the range of xbnds.}
  \item{shape}{the \emph{shape} = yheight/xwidth of the plotting regions.}
  \item{xbnds, ybnds}{horizontal and vertical limits of the binning
    region in x or y units respectively; must be numeric vector of length 2.}
  \item{xlab, ylab}{optional character strings used as labels for
    \code{x} and \code{y}.  If \code{NULL}, sensible defaults are used.}
  \item{IDs}{logical indicating if the individual cell \dQuote{IDs}
    should be returned, see also below.}
}
\value{
  an S4 object of class \code{"hexbin"}.
  It has the following slots:
  \item{cell}{vector of cell ids that can be mapped into the (x,y)
      bin centers in data units.}
  \item{count}{vector of counts in the cells.}
  \item{xcm}{The x center of mass (average of x values) for the cell.}
  \item{ycm}{The y center of mass (average of y values) for the cell.}
  \item{xbins}{ number of hexagons across the x axis. hexagon inner
      diameter =diff(xbnds)/xbins in x units}
  \item{shape}{plot shape which is yheight(inches) / xwidth(inches)}
  \item{xbnds}{x coordinate bounds for binning and plotting}
  \item{ybnds}{y coordinate bounds for binning and plotting}
  \item{dimen}{The i and j limits of cnt treated as a matrix cnt[i,j]}
  \item{n}{number of (non NA) (x,y) points, i.e., \code{sum(* @count)}.}
  \item{ncells}{number of cells, i.e., \code{length(* @count)}, etc}
  \item{call}{the function call.}
  \item{xlab, ylab}{character strings to be used as axis labels.}
  \item{cID}{of class, \code{"integer or NULL"}, only if \code{IDs}
    was true, an integer vector of length \code{n} where
    \code{cID[i]} is the cell number of the i-th original point
    \code{(x[i], y[i])}.  Consequently, the \code{cell} and \code{count}
    slots are the same as the \code{\link{names}} and entries of
    \code{table(cID)}, see the example.}
}

\seealso{
  \code{\link{hcell2xy}}%, \code{\link{hcell}},
% FIXME
  \code{\link{gplot.hexbin}},% \code{\link{hboxplot}},
%   \code{\link{hdiffplot}}, \code{\link{hmatplot}},
   \code{\link{grid.hexagons}}, \code{\link{grid.hexlegend}}.
}

\references{
  Carr, D. B. et al. (1987)
  Scatterplot Matrix Techniques for Large \eqn{N}.
  \emph{JASA} \bold{83}, 398, 424--436.
}

\details{
  Returns counts for non-empty cells only.  The plot shape must be maintained for
  hexagons to appear with equal sides.  Some calculations are in single
  precision.

  Note that when plotting a \code{hexbin} object, the
  \pkg{grid} package is used.
  You must use its graphics (or those from package \pkg{lattice} if you
  know how) to add to such plots.
}

\examples{
set.seed(101)
x <- rnorm(10000)
y <- rnorm(10000)
(bin <- hexbin(x, y))
## or
plot(hexbin(x, y + x*(x+1)/4),
     main = "(X, X(X+1)/4 + Y)  where X,Y ~ rnorm(10000)")

## Using plot method for hexbin objects:
plot(bin, style = "nested.lattice")

hbi <- hexbin(y ~ x, xbins = 80, IDs= TRUE)
str(hbi)
tI <- table(hbi@cID)
stopifnot(names(tI) == hbi@cell,
                tI  == hbi@count)

## NA's now work too:
x[runif(6, 0, length(x))] <- NA
y[runif(7, 0, length(y))] <- NA
hbN <- hexbin(x,y)
summary(hbN)
}
\keyword{dplot}
