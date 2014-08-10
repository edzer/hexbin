\name{hexGraphPaper}
\alias{hexGraphPaper}
\alias{hgridcent}
\title{Create a Hexgon Grid}
\description{
 Creates a hexagon grid that can be added to a plot created with grid
 graphics.
}
\usage{
hexGraphPaper(hb, xbnds = NULL, ybnds = NULL, xbins = 30, shape = 1,
              add = TRUE, fill.edges = 1, fill = 0, border = 1)

hgridcent(xbins, xbnds, ybnds, shape, edge.add = 0)
}
\arguments{
  \item{hb}{a object of class \code{"hexbin"}, typically produced by
    \code{\link{hexbin}(*)}.}
  \item{xbnds, ybnds}{horizontal and vertical limits of the binning
    region in x or y units respectively; must be numeric vector of
    length 2.}
  \item{xbins}{the number of bins partitioning the range of xbnds.}
  \item{shape}{the \emph{shape} = yheight/xwidth of the plotting regions.}
  \item{add}{a logical value indicating whether or not to add the grid
    to the current plot.}
  \item{fill.edges}{integer number of hexagons to add around the border}
  \item{fill}{the fill color for the hexagons}
  \item{border}{the color of the border of the hexagons}
  \item{edge.add}{offset (typically \code{fill.edges} above) used in
    \code{hgridcent}.}
}
\details{
  If a hexbin object is given then the parameters xbins and shape are
  ignored. Different bounds can still be specified. The \code{fill.edges}
  parameter should be an integer. \code{fill.edges} takes the current
  grid and adds a layer of hexagons around the grid for each level of
  fill. So for example if \code{fill.edges= 2} than the dimensions of
  the grid would be \code{(i,j)+4}.

  \code{hgridcent()} is the utility function computing the resulting
  list (see section \dQuote{Value}).

  \strong{WARNING! If using a hexVP be sure to set clip to "on", otherwise the
  hexagon grid will bleed over the plot edges.}
}
\value{
 Invisibly returns a list with th following components
  \item{x}{The x coordinates of the grid}
  \item{y}{the y coordinates of the grid}
  \item{dimen}{a vector of length 2 gining the rows and columns of the grid}
  \item{dx}{the horizontal diameter of the hexagons}
  \item{dy}{the vertical diameter of the hexagons}
}
\author{Nicholas Lewin-Koh}
\seealso{\code{\link{hcell2xy}}, \code{\link{hexpolygon}},
         \code{\link{grid.hexagons}}}
\examples{
 x <- rnorm(10000)
 y <- rnorm(10000,x,x)
 hbin <- hexbin(x,y)
 hvp <- plot(hbin,type="n")
 pushHexport(hvp$plot,clip="on")
 hexGraphPaper(hbin,border=grey(.8))
 grid.hexagons(hbin)
}
\keyword{aplot}
\keyword{dplot}
