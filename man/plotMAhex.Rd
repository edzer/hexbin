\name{plotMAhex}
\alias{plotMAhex}
\title{MA-plot using hexagon bins}
\description{
  Creates an MA-plot using hexagons with color/glyph coding for control spots.
}
\usage{
plotMAhex(MA, array = 1, xlab = "A", ylab = "M",
          main = colnames(MA)[array], xlim = NULL, ylim = NULL,
          status = NULL, values, pch, col, cex, nbin = 40,
          zero.weights = FALSE, style = "colorscale", legend = 1.2,
          lcex = 1, minarea = 0.04, maxarea = 0.8, mincnt = 2,
          maxcnt = NULL, trans = NULL, inv = NULL, colorcut = NULL,
          border = NULL, density = NULL, pen = NULL,
          colramp = function(n) { LinGray(n, beg = 90, end = 15) },
          newpage = TRUE, type = c("p", "l", "n"),
          xaxt = c("s", "n"), yaxt = c("s", "n"),
          verbose = getOption("verbose"))
}
\arguments{
  \item{MA}{an \code{RGList}, \code{MAList} or \code{MArrayLM} object,
    or any list with components \code{M} containing log-ratios and
    \code{A} containing average intensities. Alternatively a
    \code{matrix}, \code{Affybatch} or \code{ExpressionSet} object.}
  \item{array}{integer giving the array to be plotted. Corresponds to
    columns of \code{M} and \code{A}.}
  \item{xlab, ylab, main}{character strings giving label for x-axis,
    y-axis or main tile of the plot.}
  \item{xlim, ylim}{numeric vectors of length 2 giving limits for x-axis
    (or y-axis respectively), defaulting to min and max of the data.}
  \item{status}{character vector giving the control status of each spot
    on the array, of same length as the number of rows of \code{MA$M}.
    If omitted, all points are plotted in the default color, symbol and size.}
  \item{values}{character vector giving values of \code{status} to be
    highlighted on the plot. Defaults to unique values of \code{status}.
    Ignored if there is no \code{status} vector.}
  \item{pch}{vector or list of plotting characters. Default to integer code 16.
    Ignored is there is no \code{status} vector.}
  \item{col}{numeric or character vector of colors, of the same length
    as \code{values}. Defaults to \code{1:length(values)}.  Ignored if
    there is no \code{status} vector.}
  \item{cex}{numeric vector of plot symbol expansions, of the the same
    length as \code{values}.  Defaults to 0.2 for the most common status
    value and 1 for the others. Ignored if there is no \code{status} vector.}
  \item{nbin}{ Number of bins }  %% << FIXME
  \item{zero.weights}{logical, should spots with zero or negative
    weights be plotted?}
  \item{style}{string specifying the style of hexagon plot,
      see \code{\link{grid.hexagons}} for the possibilities.}
  \item{legend}{numeric width of the legend in inches of \code{FALSE}.
    In the latter case, or when \code{0}, no legend is not produced.}
  \item{lcex}{characters  expansion size for the text in the legend.}
  \item{minarea}{fraction of cell area for the lowest count.}
  \item{maxarea}{fraction of the cell area for the largest count.}
  \item{mincnt}{cells with fewer counts are ignored.}
  \item{maxcnt}{cells with more counts are ignored.}
  \item{trans}{\code{\link{function}} specifying a transformation for
    the counts such as \code{sqrt}.}
    \item{inv}{the inverse transformation of \code{trans}.}
  \item{colorcut}{vector of values covering [0, 1] that determine
    hexagon color class boundaries and hexagon legend size boundaries.
    Alternatively, an integer (\code{<= maxcnt}) specifying the
    \emph{number} of equispaced colorcut values in [0,1].}
  \item{border, density, pen}{color for polygon borders and filling of
    each hexagon drawn, passed to \code{\link{grid.hexagons}}.}
  \item{colramp}{function accepting an integer \code{n} as an argument and
    returning n colors.}
  \item{newpage}{should a new page start?}
  \item{type, xaxt, yaxt}{strings to be used (when set to \code{"n"}) for
    suppressing the plotting of hexagon symbols, or the x- or y-axis,
    respectively.}
  \item{verbose}{logical indicating if some diagnostic output should happen.}
}

\details{
  An MA-plot is a plot of log-intensity ratios (M-values) versus
  log-intensity averages (A-values). If \code{MA} is an \code{RGList} or
  \code{MAList} then this function produces an ordinary within-array
  MA-plot. If \code{MA} is an \code{MArrayLM} object, then the plot is an
  fitted model MA-plot in which the estimated coefficient is on the y-axis
  and the average A-value is on the x-axis.

  If \code{MA} is a \code{matrix} or \code{ExpressionSet} object, then this
  function produces a between-array MA-plot. In this case the A-values in
  the plot are the average log-intensities across the arrays and the
  M-values are the deviations of the log-intensities for the specified
  array from the average. If there are more than five arrays, then the
  average is computed robustly using medians. With five or fewer arrays,
  it is computed by means.

  The \code{status} vector is intended to specify the control status of
  each spot, for example "gene", "ratio control", "house keeping gene",
  "buffer" and so on. The vector is usually computed using the function
  \code{\link[limma]{controlStatus}} from package \pkg{limma} and a
  spot-types file.  However the function may be used to highlight any
  subset of spots.

  The arguments \code{values}, \code{pch}, \code{col} and \code{cex}
  can be included as attributes to \code{status} instead of being
  passed as arguments to \code{plotMA}.

  See \code{\link[graphics]{points}} for possible values for \code{pch},
  \code{col} and \code{cex}.
}

\value{
  A plot is created on the current graphics device.
  and a list with the following items is returned invisibly:
  \item{plot.vp}{the \code{\link{hexViewport}} constructed and used.}
  \item{legend.vp}{if a legend has been produced, its
    \code{\link[grid]{viewport}}.}
  \item{hbin}{a \code{hexbin} object built with A as the x coordinate
    and M as the y coordinate.}
}

\references{See \url{http://www.statsci.org/micrarra/refs/maplots.html}}

\author{Nicholas Lewin-Koh, adapted from code by Gordon Smyth}

\seealso{\code{\link[limma:plotma]{plotMA}} from package \pkg{limma},
  and \code{\link{gplot.hexbin}}.
}

\examples{
  if(require(marray)){ %% for the data only --> data(swirl, package="marray")
    data(swirl)
    hb <- plotMAhex(swirl[,1],newpage=FALSE,
                    main = "M vs A plot with hexagons", legend=0)
    hexVP.abline(hb$plot.vp,h=0,col=gray(.6))
    hexMA.loess(hb)
  }
}
\keyword{hplot}
