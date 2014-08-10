\name{hboxplot}
\alias{hboxplot}
\title{2-D Generalization of Boxplot}

\description{
  If \code{bin} is an \emph{eroded} \code{\link{hexbin}} object, i.e.,
  an \code{erodebin} object, \code{hboxplot()} plots the high counts cells
  selected by \code{\link{erode}()}.  By default, the high counts
  cells contain 50 percent of the counts so analagous to the
  interquartile \dQuote{range}.  The function distinguishes the last
  cells eroded using color.  These cells correspond to one definition of the
  bivariate median.
%% FIXME ^^ (bad style, content +- ok)
}

\usage{
hboxplot(bin, xbnds = NULL, ybnds = NULL,
         density, border = c(0, grey(0.7)), pen = c(2, 3),
         unzoom = 1.1, clip ="off", reshape = FALSE,
         xlab = NULL, ylab = NULL, main = "")
}
\arguments{
  \item{bin}{an object of class \code{\link{hexbin}}.}
  \item{xbnds,ybnds}{global x- and y-axis plotting limits for multiple
    plots.}
  \item{density, border}{arguments for \code{\link{polygon}()} each of
    length two, the first for the median, the second for the other cells.}
  \item{pen}{colors (\dQuote{pen numbers}) for \code{polygon()}.}
  \item{unzoom}{plot limit expansion factor when \code{xbnds} is
    missing.}
  \item{clip}{either 'on' or 'off' are the allowed arguments, when on
    everything is clipped to the plotting region.}
  \item{reshape}{logical value to reshape the plot although \code{xbnds}
    and \code{ybnds} are present.}
  \item{xlab, ylab, main}{x- and y- axis labels and main title}
}

\value{
  invisibly, the \code{\link{hexViewport}()} used internally.
  Used to add to the plot afterwards.
}

\references{ see in \code{\link{grid.hexagons}}.}

\details{
  The \code{density}, \code{border}, and \code{pen} arguments correspond
  to the \code{\link{polygon}} function calls for plotting two types of
  cells.  The cell types, pen numbers and suggested colors are\cr
  \tabular{lll}{
    TYPE				   \tab PEN \tab COLOR \cr
    cells of bin                           \tab 2  \tab light gray \cr
    last eroded cells of bin (median cells)\tab 1  \tab black \cr
  }

  The erode components of the hexbin objects must be present for the
  medians cells to plot.

  When \code{xbnds} is missing or \code{reshape} is true, the plot
  changes graphics parameters and resets them.   When \code{xbnds} is
  missing the function also zooms in based on the available data to
  provide increased resolution.

  The zoom used the hexagon cell centers.  The unzoom argument backs off
  a bit so the whole hexagon will fit in the plot.

  \code{Hboxplot()} is used as a stand alone function, for producing separate
  legends .....%%FIXME for \code{\link{hmatplot}()} and for panels in
  %% \code{\link{hmatplot}()}.
}

\seealso{
  \code{\link{hexbin}}, \code{\link{erode}},
  %\code{\link{smooth.hexbin}},
  \code{\link{hcell2xy}},%  \code{\link{hcell}},
  \code{\link{gplot.hexbin}},
%  \code{\link{hmatplot}},
  \code{\link{grid.hexagons}}, \code{\link{grid.hexlegend}}
}

\examples{
\dontshow{set.seed(753)}
##  boxplot of smoothed counts
x <- rnorm(10000)
y <- rnorm(10000)

bin <- hexbin(x,y)
erodebin <- erode(smooth.hexbin(bin))

hboxplot(erodebin)
hboxplot(erodebin, density = c(32,7), border = c(2,4))
hp <- hboxplot(erodebin, density = c(-1,17),
               main = "hboxplot(erode*(smooth*(.)))")
pushHexport(hp)
library("grid")
grid.points(x[1:10], y[1:10])# just non-sense to show the principle
popViewport()
}
\keyword{hplot}
