\name{hexVP-class}
\docType{class}
\alias{hexVP-class}
\alias{getFig,hexVP-method}
\alias{getMargins,hexVP-method}
\alias{getPlt,hexVP-method}
\alias{getXscale,hexVP-method}
\alias{getYscale,hexVP-method}

\title{Formal class "hexVP" of a Hexagon Viewport}
\description{
  Hexagon Viewports are \dQuote{value-added} grid viewports (see
  \code{\link[grid]{viewport}}) where the extra slots contain scaling and
  \dQuote{embedding} information. A hexViewport is created my taking the
  available area in the cuurent viewport on the graphics device and
  maximizing the amount of area with a fied aspect ratio. The default
  when the shape parameter is 1, is a 1:1 aspect ratio in terms of the
  size of the viewport, not the scale of the x and y axis. The plotting
  area is centered within the existing margins and the maximum size
  determined. Extra area is then allocated to the margins. This viewport
  is replicated twice, once with clipping set to "on" and once with
  clipping "off". This feature can be used for toggling clipping on and
  off while editing the plot.
}
\section{Objects from the Class}{
  Objects are typically created by calls to \code{\link{hexViewport}()}
  or by low level calls of the form \code{new("hexVP", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{hexVp.off}:}{Object of class \code{"viewport"} with
      clipping set to off, see \code{\link[grid]{viewport}}.}
    \item{\code{hexVp.on}:}{Object of class \code{"viewport"}, with the same
      dimensions and parameters as hexVp.off, but with
      clipping set to on, see \code{\link[grid]{viewport}}.}
    \item{\code{hp.name}:}{The name of the viewport for searching a vptree}.
    \item{\code{mar}:}{\code{\link[grid]{unit}} vector of four margins
      (typically in \code{"lines"}).}
    \item{\code{fig}:}{\code{\link[grid]{unit}} vector of two figure sizes
      (typically in \code{"npc"}).}
    \item{\code{plt}:}{\code{\link[grid]{unit}} vector of two figure sizes
      (typically in \code{"npc"}).}
    %% MM {FIXME?}: Is n't this simply ``xlim'' - then call it so!
    %% NL, yes it is, but xscale and yscale is the parameters used by
    %% grid.
    \item{\code{shape}:}{The shape parameter from the plotted
      \code{\link[hexbin]{hexbin}} object.}
    \item{\code{xscale}:}{numeric of length two specifying x-range.}
    \item{\code{yscale}:}{numeric of length two specifying y-range.}
  }
}
\section{Methods}{
  These are methods accessing the slots of corresponding name.
  \describe{
    \item{getFig}{\code{signature(hvp = "hexVP")}: ... }
    \item{getMargins}{\code{signature(hvp = "hexVP")}: ... }
    \item{getPlt}{\code{signature(hvp = "hexVP")}: ... }
    \item{getXscale}{\code{signature(hvp = "hexVP")}: ... }
    \item{getYscale}{\code{signature(hvp = "hexVP")}: ... }
  }
}
\author{
  Nicholas Lewin-Koh \email{kohnicho@comp.nus.edu.sg}.
}
\seealso{
  The constructor function \code{\link{hexViewport}}.
  \code{\link{hexbin}}, and its S4 plotting method,
  \code{\link{gplot.hexbin}}.
}
\examples{
 library("grid")
 example(hexViewport, echo=FALSE)
 ## continued:
 str(P$plot.vp)
}
\keyword{classes}
