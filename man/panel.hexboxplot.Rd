\name{panel.hexboxplot}
\alias{panel.hexboxplot}
\title{Boxplot for hexbin lattice plot}
\description{
A panel function to add a boxplot to a hexbin lattice plot.
}
\usage{
panel.hexboxplot(x, y, xbins = 30,
                 xbnds = c("data", "panel"), ybnds = c("data", "panel"),
                 .prelim = FALSE, .cpl = current.panel.limits(),
                 .xlim = .cpl$xlim, .ylim = .cpl$ylim,
                 .aspect.ratio, type = character(0), cdfcut = 0.25,
                 shadow = 0.05, ..., check.erosion = TRUE)
}
\arguments{
  \item{x, y}{numeric vector or factor.}
  \item{xbins}{the number of bins partitioning the range of xbnds.}
  \item{xbnds, ybnds}{horizontal and vertical limits of the binning
    region in x or y units respectively; must be numeric vector of
    length 2.}
  \item{.prelim, .cpl, .xlim, .ylim, .aspect.ratio}{for internal use.}
  \item{type}{character vector controlling additional augmentation of
    the display. A \code{"g"} in \code{type} adds a reference grid, an
    \code{"hg"} adds a hexagonal grid.}
  \item{cdfcut}{number in (0,1) indicating the confidence level for the
    erosion limits.  See \code{\link{erode.hexbin}} for more information.}
  \item{shadow}{number in (0,1) indicating the confidence level for the
    erosion limits of a boxplot shadow.  See \code{\link{erode.hexbin}}
    for more information.}
  \item{\dots}{potential further arguments passed on.}
  \item{check.erosion}{logical indicating only eroded points should be
    used for \code{"erodebin"} objects; simply passed to
    \code{\link{hcell2xy}}, see its documentation.}
}
\value{
 There is no return value from this function. The results are plotted on
 the current active device.
}
\author{Nicholas Lewin-Koh \email{nikko@hailmail.net}}
\seealso{\code{\link{hexbinplot}}, \code{\link{panel.hexgrid}},
  \code{\link[lattice]{panel.bwplot}}
}
\examples{
mixdata <-
    data.frame(x = c(rnorm(5000),rnorm(5000,4,1.5)),
               y = rep(1:2, 5000))
hexbinplot(y ~ x, mixdata, panel = panel.hexboxplot)
}
\keyword{hplot}
