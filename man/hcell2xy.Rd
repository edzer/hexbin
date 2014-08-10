\name{hcell2xy}
\alias{hcell2xy}
\alias{hcell2xy,hexbin-method}

\title{Compute X and Y Coordinates for Hexagon Cells}

\description{
  Computes x and y coordinates from hexagon cell id's.
}
\usage{
hcell2xy(hbin, check.erosion = TRUE)
}
\arguments{
  \item{hbin}{a object of class \code{"hexbin"}, typically produced by
    \code{\link{hexbin}(*)}.}
  \item{check.erosion}{logical indicating if only the eroded points
    should be returned in the case where \code{hbin} inherits from
    \code{"erodebin"} (see \code{\link{erodebin-class}}); is \code{TRUE}
    by default.}
}
\value{
  A list with two components of the same length as \code{bin$cell},
  \item{x}{}
  \item{y}{}
}

%%FIXME \references{see in \code{\link{hcell}}.}

\details{
  The hexbin object \code{hbin} contains all the needed information.
  The purpose of this function is to reduce storage.  The cost is
  additional calculation.
}

\seealso{%%FIXME \code{\link{hcell}}, \code{\link{hray}},
  \code{\link{hexbin}}.
}

\examples{
x <- rnorm(10000)
y <- rnorm(10000)
plot(x,y, pch=".")
hbin <- hexbin(x,y)
str(xys <- hcell2xy(hbin))
points(xys, cex=1.5, col=2) ; title("hcell2xy( hexbin(..) )", col.main=2)

%% __________ FIXME ________
\dontshow{
## Temporal trends with confidence bounds plotted on a map:
## <NOT YET> Illustration only pending access to user functions
##      mtapply()            # like tapply but for matrices
##      sens.season.slope()  # computes sen's seasonal slope

## This part does not work and commented out
#hbin  <- hcell(dat$x,dat$y)   #  x and y are in map projection units
#newdat < dat[,c('month','year','value')]    # extract columns
#stats <- mtapply(newdat,bin$cell,sens.season.slope,season=12)
#plot(mymap,type='l')       # map boundaries in map projection units
#xy <- hcell2xy(hbin)  # x and y coordinates for hexagon cell centers
#hray(xy$x, xy$y,val=stat[,1],lower= stat[,2],upper=stat[,3])
}
}
\keyword{manip}
