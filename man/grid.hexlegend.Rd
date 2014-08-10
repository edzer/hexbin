\name{grid.hexlegend}
\alias{grid.hexlegend}
\title{Add a Legend to a Hexbin Plot}
\description{
  Plots the legend for the \code{plot} method of \code{\link{hexbin}}.
  Provides a legend indicating the count representations.
}
\usage{
grid.hexlegend(legend, ysize, lcex, inner, style = ,
  	   minarea = 0.05, maxarea = 0.8, mincnt = 1, maxcnt, trans = NULL,
	   inv = NULL, colorcut, density = NULL, border = NULL, pen = NULL,
	   colramp = function(n) { LinGray(n,beg = 90,end = 15) },
           leg.unit = "native")
}
\arguments{
  \item{legend}{positive number giving width of the legend in inches.}
  \item{ysize}{height of legend in inches}
  \item{lcex}{the characters expansion size for the text in the legend,
    see \code{\link{par}(cex=)}.}
  \item{inner}{the inner diameter of a hexagon in inches.}
  \item{style}{the hexagon style; see \code{\link{grid.hexagons}}.}
  \item{minarea, maxarea}{fraction of the cell area for the lowest and largest
    count, respectively.}
  \item{mincnt, maxcnt}{minimum and maximum count accepted in \code{plot}.}
  \item{trans}{a transformation function for the counts such as
    \code{\link{sqrt}}.}
  \item{inv}{the inverse transformation function.}
  \item{colorcut}{numeric vector of values covering [0, 1] the determine
    hexagon color classes boundaries and hexagon legend size boundaries.}
  \item{border}{argument for \code{\link{polygon}()}.  Draw the border
    for each hexagon.}
  \item{density}{argument for \code{\link{polygon}()} filling.  A
    \code{0} causes the polygon not to be filled.}
  \item{pen}{color argument used for \code{\link{polygon}(col = .)}.
    Determines the color with which the polygon will be filled.}
  \item{colramp}{function accepting an integer \code{n} as an argument and
    returning n colors.}
  \item{leg.unit}{unit to use}%FIXME

}
\details{
  The \code{plot} method for \code{\link{hexbin}} objects calls this function
  to produce a legend
  by setting the graphics parameters, so \code{hex.legend} itself is not a
  standalone function.

  The legend function is \bold{preliminary}.  Later version will include
  refinements and handle extreme cases (small and large) for cell size
  and counts.

  See the \bold{Details} section of \code{\link{grid.hexagons}}'s help page.
}
\value{
  This function does not return any value.
}
\references{ see in \code{\link{grid.hexagons}}.}
\author{
  Dan Carr <dcarr@voxel.galaxy.gmu.edu>

  ported by Nicholas Lewin-Koh <kohnicho@comp.nus.edu.sg>
}

\seealso{\code{\link{hexbin}}, \code{\link{grid.hexagons}},
% FIXME
   \code{\link{smooth.hexbin}}, \code{\link{erode.hexbin}},
%   \code{\link{hcell}},
  \code{\link{hcell2xy}},
   \code{\link{gplot.hexbin}},% \code{\link{hboxplot}},% \code{\link{hdiffplot}},
%   \code{\link{hmatplot}}
}

\examples{
## Not a stand alone function; typically only called from plot.hexbin()
%% Hence we should not run it here!
%% FIXME: Improve hex.legend() such that it *can* be added to plots !!!
\dontrun{
  grid.hexlegend(legend = 2, ysize = 1,lcex=8,inner=0.2,
                 maxcnt = 100, colorcut = c(0.5,0.5))
}
}
\keyword{aplot}
