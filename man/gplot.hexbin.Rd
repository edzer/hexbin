\name{gplot.hexbin}
\alias{gplot.hexbin}
\alias{plot,hexbin,missing-method}
\title{Plotting Hexagon Cells with a Legend}
\description{
  Plots Hexagons visualizing the counts in an hexbin object.  Different
  styles are availables.  Provides a legend indicating the count
  representations.
}
\usage{
%% In future: No longer export   gplot.hexbin() !
gplot.hexbin(x, style = "colorscale", legend = 1.2, lcex = 1,
      minarea = 0.04, maxarea = 0.8, mincnt = 1, maxcnt = max(x@count),
      trans = NULL, inv = NULL, colorcut = seq(0, 1, length = min(17, maxcnt)),
      border = NULL, density = NULL, pen = NULL,
      colramp = function(n) LinGray(n,beg = 90,end = 15),
      xlab = "", ylab = "", main = "", newpage = TRUE,
      type = c("p", "l", "n"), xaxt = c("s", "n"), yaxt = c("s", "n"),
      clip = "on", verbose = getOption("verbose"))
%% FIXME: This is the S4 plot method for 'hexbin'
%%        currently also exported "standalone" - for testing,debugging..
%% we'd really don't want to repeat the argument list; use \synopsis{.} ?
\S4method{plot}{hexbin,missing}(x, style = "colorscale", legend = 1.2, lcex = 1,
      minarea = 0.04, maxarea = 0.8, mincnt = 1, maxcnt = max(x@count),
      trans = NULL, inv = NULL, colorcut = seq(0, 1, length = min(17, maxcnt)),
      border = NULL, density = NULL, pen = NULL,
      colramp = function(n) LinGray(n,beg = 90,end = 15),
      xlab = "", ylab = "", main = "", newpage = TRUE,
      type = c("p", "l", "n"), xaxt = c("s", "n"), yaxt = c("s", "n"),
      clip = "on", verbose = getOption("verbose"))
}
\arguments{
  \item{x}{an object of class \code{\link{hexbin}}.}
%   \item{y}{(required by the S4 method for \code{\link{plot}} but unused
%     here; must be missing)}
  \item{style}{string specifying the style of hexagon plot,
      see \code{\link{grid.hexagons}} for the possibilities.}
  \item{legend}{numeric width of the legend in inches of \code{FALSE}.
    In the latter case, or when \code{0}, no legend is not produced.}
  \item{lcex}{characters  expansion size for the text in the legend}
  \item{minarea}{fraction of cell area for the lowest count}
  \item{maxarea}{fraction of the cell area for the largest count}
  \item{mincnt}{cells with fewer counts are ignored.}
  \item{maxcnt}{cells with more counts are ignored.}
  \item{trans}{\code{\link{function}} specifying a transformation for
    the counts such as \code{sqrt}.}
% FIXME: use better description of these in hexagons()  -- or use same
% ----   help page ?!
  \item{inv}{the inverse transformation of \code{trans}.}
  \item{colorcut}{vector of values covering [0, 1] that determine
    hexagon color class boundaries and hexagon legend size boundaries.
    Alternatively, an integer (\code{<= maxcnt}) specifying the
    \emph{number} of equispaced colorcut values in [0,1].}
  \item{border, density, pen}{color for polygon borders and filling of
    each hexagon drawn, passed to \code{\link{grid.hexagons}}.}
  \item{colramp}{function accepting an integer \code{n} as an argument and
    returning n colors.}
  \item{xlab, ylab}{x- and y-axis label.}
  \item{main}{main title.}
  \item{newpage}{should a new page start?.}
  \item{type, xaxt, yaxt}{strings to be used (when set to \code{"n"}) for
    suppressing the plotting of hexagon symbols, or the x- or y-axis,
    respectively.}
  \item{clip}{either 'on' or 'off' are the allowed arguments, when on
    everything is clipped to the plotting region.}
  \item{verbose}{logical indicating if some diagnostic output should happen.}
  \item{xlim, ylim}{numeric vectors of length 2, giving the x and y coordinates ranges.}
}
\details{
  This is the (S4) \code{\link{plot}} method for \code{\link{hexbin}} (and
  \code{erodebin}) objects (\link{erodebin-class}).

  To use the standalone function
  \code{gplot.hexbin()} is \bold{\emph{deprecated}}.
  For \code{style}, \code{minarea} etc, see the \bold{Details} section of
  \code{\link{grid.hexagons}}'s help page.

  The legend functionality is somewhat preliminary.  Later versions may
  include refinements and handle extreme cases (small and large) for
  cell size and counts.

  All arguments of \code{gplot.hexbin} can also be used for the S4 
  \code{\link{plot}} method.
}
\value{
  invisibly, a list with components
  \item{plot.vp}{the \code{\link{hexViewport}} constructed and used.}
  \item{legend.vp}{if a legend has been produced, its
    \code{\link[grid]{viewport}}.}
}

\references{ see in \code{\link{grid.hexagons}}.}
\author{
  Dan Carr \email{dcarr@voxel.galaxy.gmu.edu},
  ported by Nicholas Lewin-Koh \email{kohnicho@comp.nus.edu.sg} and
  Martin Maechler.
}
\seealso{\code{\link{hexbin}}, \code{\link{hexViewport}},
  \code{\link{smooth.hexbin}},
  \code{\link{erode.hexbin}},
  \code{\link{hcell2xy}}, \code{\link{hboxplot}},
  \code{\link{hdiffplot}}.  %%, \code{\link{hmatplot}}.
}
\examples{
## 1) simple binning of spherical normal:
x <- rnorm(10000)
y <- rnorm(10000)
bin <- hexbin(x,y)

## Plot method for hexbin !
## ---- ------     --------
plot(bin)
# nested lattice
plot(bin, style= "nested.lattice")

# controlling the colorscheme
plot(bin, colramp=BTY, colorcut=c(0,.1,.2,.3,.4,.6,1))

## 2) A mixture distribution
x <- c(rnorm(5000),rnorm(5000,4,1.5))
y <- c(rnorm(5000),rnorm(5000,2,3))
bin <- hexbin(x,y)

pens <- cbind(c("#ECE2F0","#A6BDDB","#1C9099"),
              c("#FFF7BC","#FEC44F","#D95F0E"))
plot(bin, style = "nested.lattice", pen=pens)
# now really crazy
plot(bin, style = "nested.lattice", pen=pens,border=2,density=35)

# lower resolution binning and overplotting with counts
bin <- hexbin(x,y,xbins=25)
P <- plot(bin, style="lattice", legend=FALSE,
          minarea=1, maxarea=1, border="white")
##
%% FIXME!
library("grid")
pushHexport(P$plot.vp)
xy <- hcell2xy(bin)
  # to show points rather than counts :
grid.points(x,y,pch=18,gp=gpar(cex=.3,col="green"))
grid.text(as.character(bin@count), xy$x,xy$y,
          gp=gpar(cex=0.3, col="red"),default.units="native")
popViewport()

# Be creative, have fun!
}
\keyword{hplot}
