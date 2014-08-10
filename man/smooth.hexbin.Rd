\name{smooth.hexbin}
\alias{smooth.hexbin}
\alias{smoothbin-class}
\title{Hexagon Bin Smoothing}
\description{
  Given a \code{"hexbin"} (hexagon bin) object, compute a discrete
  kernel smoother that covers seven cells, namely a center cell and its
  six neighbors.  With two iterations the kernel effectively covers
  1+6+12=19 cells.
}
\usage{
smooth.hexbin(bin, wts=c(48,4,1))
}
\arguments{
  \item{bin}{object of class \code{"hexbin"}, typically resulting from
    \code{\link{hexbin}()} or \code{\link{erode,hexbin-method}}.}
  \item{wts}{numeric vector of length 3 for relative weights of the
    center, the six neighbor cells, and twelve second neighbors.}
}
\value{
  an object of class \code{"smoothbin"}, extending class
  \code{"hexbin"}, see \code{\link{hexbin}}.
  The object includes the additional slot \code{wts}.
}
\references{see \code{\link{grid.hexagons}} and \code{\link{hexbin}}.}

\details{
  This discrete kernel smoother uses the center cell, immediate
  neighbors and second neighbors to smooth the counts.  The counts for
  each resulting cell is a linear combination of previous cell counts
  and weights. The weights are
  \tabular{ll}{
    1 center cell,       \tab weight = wts[1]\cr
    6 immediate neighbors\tab weight = wts[2]\cr
    12 second neighbors  \tab  weight =wts[3]\cr
  }
  If a cell, its immediate and second neighbors all have a value of
  \code{max(cnt)}, the new maximum count would be
  \code{max(cnt)*sum(wts)}.  It is possible for the counts to overflow.

  The domain for cells with positive counts increases.  The hexbin
  slots  \code{xbins}, \code{xbnds}, \code{ybnds}, and \code{dimen} all
  reflect this increase.
  Note that usually \code{dimen[2] = xbins+1}.

  The intent was to provide a fast, iterated, immediate neighbor smoother.
  However, the current hexbin plotting routines only support shifting 
  even numbered rows to the right.  Future work can

  (1) add a shift indicator to hexbin objects that indicates left or
      right shifting.\cr
  (2) generalize plot.hexbin() and hexagons()\cr
  (3) provide an iterated kernel.\cr

  With \code{wts[3]=0}, the smoother only uses the immediate neighbors.
  With a shift indicator the domain could increase by 2 rows (one bottom
  and on top) and 2 columns (one left and one right).  However the current
  implementation increases the domain by 4 rows and 4 columns, thus
  reducing plotting resolution.
}

\seealso{
  \code{\link{hexbin}}, \code{\link{erode.hexbin}}, %MISSING \code{\link{hthin}},
  \code{\link{hcell2xy}},% \code{\link{hcell}},
  \code{\link{gplot.hexbin}}, \code{\link{hboxplot}},
  %\code{\link{hdiffplot}}, \code{\link{hmatplot}},
  \code{\link{grid.hexagons}}, \code{\link{grid.hexlegend}}.
}
\examples{
x <- rnorm(10000)
y <- rnorm(10000)
bin <- hexbin(x,y)
# show the smooth counts in gray level
smbin  <- smooth.hexbin(bin)
plot(smbin, main = "smooth.hexbin(.)")

# Compare the smooth and the origin
smbin1 <- smbin
smbin1@count <- as.integer(ceiling(smbin@count/sum(smbin@wts)))
plot(smbin1)
smbin2 <- smooth.hexbin(bin,wts=c(1,0,0))  # expand the domain for comparability
plot(smbin2)
}
\keyword{misc}
