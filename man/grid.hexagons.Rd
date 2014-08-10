\name{grid.hexagons}
\alias{grid.hexagons}
\title{Add Hexagon Cells to Plot}
\description{
  Plots cells in an hexbin object.  The function distinquishes among
  counts using 5 different styles.  This function is the hexagon
  plotting engine from the \code{plot} method for \code{\link{hexbin}}
  objects.
}
\usage{
grid.hexagons(dat, style = c("colorscale", "centroids", "lattice",
	      "nested.lattice", "nested.centroids", "constant.col"),
         use.count=TRUE, cell.at=NULL,
	 minarea = 0.05, maxarea = 0.8, check.erosion = TRUE,
	 mincnt = 1, maxcnt = max(dat@count), trans = NULL,
	 colorcut = seq(0, 1, length = 17),
	 density = NULL, border = NULL, pen = NULL,
	 colramp = function(n){ LinGray(n,beg = 90, end = 15) },
	 def.unit=  "native",
	 verbose = getOption("verbose"))
}
\arguments{
  \item{dat}{an object of class \code{hexbin}, see \code{\link{hexbin}}.}
  \item{style}{character string specifying the type of plotting; must be (a
    unique abbrevation) of the values given in \sQuote{Usage} above.}
  \item{use.count}{logical specifying if counts should be used.}
  \item{cell.at}{numeric vector to be plotted instead of counts, must
    besame length as the number of cells.}
  \item{minarea}{numeric, the fraction of cell area for the lowest count.}
  \item{maxarea}{the fraction of the cell area for the largest count.}
  \item{check.erosion}{logical indicating only eroded points should be
    used for \code{"erodebin"} objects; simply passed to
    \code{\link{hcell2xy}}, see its documentation.}
  \item{mincnt}{numeric; cells with counts smaller than \code{mincnt}
    are not shown.}
  \item{maxcnt}{cells with counts larger than this are not shown.}
  \item{trans}{a transformation function (or \code{NULL}) for the counts,
    e.g., \code{\link{sqrt}}.}
  \item{colorcut}{a vector of values covering [0, 1] which determine
    hexagon color class boundaries or hexagon size boundaries -- for
    \code{style = "colorscale"} only.}
  \item{density}{\code{\link[grid]{grid.polygon}} argument for shading.  0 causes
    the polygon not to be filled. \emph{This is not implemented (for
      \code{\link[grid]{grid.polygon}}) yet}.}
  \item{border}{\code{\link[grid]{grid.polygon}()} argument.  Draw the border for
    each hexagon.}
  \item{pen}{colors for \code{\link[grid]{grid.polygon}()}.  Determines the color
    with which the polygon will be filled.}
  \item{colramp}{function of an integer argument \code{n} returning n
    colors. \code{n} is determined }%% how?  FIXME
  \item{def.unit}{default \code{\link[grid]{unit}} to be used.}% FIXME
  \item{verbose}{logical indicating if some diagnostic output should happen.}
}
\section{SIDE EFFECTS}{Adds hexagons to the plot.}

\details{
  The six plotting styles have the following effect:
  \describe{
    \item{\code{style="lattice"} or \code{"centroids"}:}{

      Plots the hexagons in different sizes based on counts.  The
      \code{"lattice"} version centers the hexagons at the cell centers
      whereas \code{"centroids"} moves the hexagon centers close to the
      center of mass for the cells.  In all cases the hexagons will not
      plot outside the cell unless \code{maxarea > 1}.  Counts are rescaled
      into the interval [0,1] and colorcuts determine the class
      boundaries for sizes and counts. The pen argument for this style
      should be a single color or a vector of colors of
      \code{length(bin@count)}.}

    \item{\code{style="colorscale"}:}{
      Counts are rescaled into the interval [0,1] and colorcuts determines
      the class boundaries for the color classes.  For this style, the
      function passed as \code{colramp} is used to define the n colors for
      the n+1 color cuts. The pen argument is ignored.
      %% S-plus: In motif color options try polygon:  black 16 white
      See \code{\link{LinGray}} for the default \code{colramp} and
      alternative \dQuote{color ramp} functions.
    }
    \item{\code{style="constant.col"}:}{
      This is an even simpler alternative to \code{"colorscale"},
      using constant colors (determined \code{pen} optionally).
    }

    \item{\code{style="nested.lattice"} and \code{"nested.centroids"}:}{
      Counts are partitioned into classes by power of 10.  The encoding
      nests hexagon size within powers of 10 color contours.

      If the pen argument is used it should be a matrix of colors with 2
      columns and either \code{ceiling(log10(max(bin@count)))} or
      \code{length(bin@count)} rows.  The default uses the \R color palatte
      so that pens numbers 2-11 determine colors for completely filled
      cell Pen 2 is the color for 1's, Pen 3 is the color for 10's, etc.
      Pens numbers 12-21 determine the color of the foreground hexagons. The
      hexagon size shows the relative count for the power of 10. Different
      color schemes give different effects including 3-D illusions
      %% S-plus :
      %%   One motif color option for the first 4 powers is black \#BBB \#36F
      %%   \#0E3 \#F206 \#FFF4 \#FFF
      %%
      %%   A second option is for the first 5 power is black \#FFF \#08F \#192
      %%   \#F11 \#FF04 \#000 \#999 \#5CF \#AFA \#FAAF \#000
    }
  }

  \emph{Hexagon size encoding \code{minarea} and \code{maxarea}}
  determine the area of the smallest and largest hexagons
  plotted.  Both are expressed fractions of the bin cell size.  Typical
  values might be .04 and 1.  When both values are 1, all plotted
  hexagons are bin cell size, if \code{maxarea} is greater than 1 than
  hexagons will overlap. This is sometimes interesting with the lattice
  and centroid styles.

  \emph{Count scaling}

  \code{relcnt <- (trans(cnt)-trans(mincnt)) / (trans(maxcnt)-trans(mincnt))}
  \cr
  \code{area <- minarea + relcnt*maxarea}

  By default the transformation \code{trans()} is the identity
  function.  The legend routine requires the transformation inverse
  for some options.

  \emph{Count windowing \code{mincnt} and \code{maxcnt}}
  Only routine only plots cells with cnts in [mincnts,   maxcnts]
}
\references{
  Carr, D. B. (1991)
  Looking at Large Data Sets Using Binned Data Plots,
  pp. 7--39 in \emph{Computing and Graphics in Statistics};
  Eds. A. Buja and P. Tukey, Springer-Verlag, New York.
}
\author{
  Dan Carr <dcarr@voxel.galaxy.gmu.edu>;
  ported and extended by Nicholas Lewin-Koh \email{nikko@hailmail.net}.
}
\seealso{\code{\link{hexbin}}, \code{\link{smooth.hexbin}},
  \code{\link{erode.hexbin}}, \code{\link{hcell2xy}},% \code{\link{hcell}},
  \code{\link{gplot.hexbin}}, \code{\link{hboxplot}}, \code{\link{hdiffplot}},
  \code{\link{grid.hexlegend}}%  \code{\link{hmatplot}}
}

\examples{
library("grid")
set.seed(506)
x <- rnorm(10000)
y <- rnorm(10000)

# bin the points
bin <- hexbin(x,y)

# Typical approach uses plot( <hexbin> ) which controls the plot shape :
plot(bin, main = "Bivariate rnorm(10000)")

## but we can have more manual control:

# A mixture distribution
x <- c(rnorm(5000),rnorm(5000,4,1.5))
y <- c(rnorm(5000),rnorm(5000,2,3))
hb2 <- hexbin(x,y)

# Show color control and overplotting of hexagons
## 1) setup coordinate system:
P <- plot(hb2, type="n", main = "Bivariate mixture (10000)")# asp=1

## 2) add hexagons (in the proper viewport):
pushHexport(P$plot.vp)
grid.hexagons(hb2, style= "lattice", border = gray(.1), pen = gray(.6),
              minarea = .1, maxarea = 1.5)
popViewport()

## How to treat 'singletons' specially:
P <- plot(hb2, type="n", main = "Bivariate mixture (10000)")# asp=1
pushHexport(P$plot.vp)
grid.hexagons(hb2, style= "nested.centroids", mincnt = 2)# not the single ones
grid.hexagons(hb2, style= "centroids", maxcnt = 1, maxarea=0.04)# single points
popViewport()


%% FIXME --- this would mix grid- and traditional-graphics
%% ----- would need grid-graphics for 'gpclib' -- aaargs...
% # And if we had all the information...
% if(require(gpclib)){
%   h1 <- chull(x[1:5000], y[1:5000])
%   h2 <- chull(x[5001:10000], y[5001:10000])
%   h2 <- h2+5000
%   h1 <- as(cbind(x[1:5000],y [1:5000])[h1, ], "gpc.poly")
%   h2 <- as(cbind(x,y)[h2, ], "gpc.poly")
%   plot(hb2, type="n", main = "Bivariate mixture (10000)")# asp=1
%
%   plot(h1,poly.args = list(col ="#CCEBC5"),add = TRUE)
%   plot(h2,poly.args = list(col ="#FBB4AE"),add = TRUE)
%   plot(intersect(h1, h2), poly.args = list(col = 2), add = TRUE)
%   grid.hexagons(hb2, style= "centroids", border = gray(.1), pen = gray(.6),
%                 minarea = .1, maxarea = 1.5)
% }

}
\keyword{aplot}
