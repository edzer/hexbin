\name{hexplom}
\title{Hexbin Plot Matrices}
\alias{hexplom}
\alias{hexplom.formula}
\alias{hexplom.data.frame}
\alias{hexplom.matrix}
\alias{panel.hexplom}

\usage{
hexplom(x, data, \dots)

\method{hexplom}{formula}(x, data = NULL, \dots)

\method{hexplom}{data.frame}(x, data = NULL, \dots, groups = NULL,
   subset = TRUE) 

\method{hexplom}{matrix}(x, data = NULL, \dots, groups = NULL, subset = TRUE)

panel.hexplom(\dots)
}

\description{
  \code{hexplom} draws Conditional Hexbin Plot Matrices.  It is similar
  to \code{splom}, expect that the default display is different.
  Specifically, the default display is created using
  \code{panel.hexplom}, which is an alias for \code{panel.hexbinplot}.
}

\arguments{
  \item{x}{
    The object on which method dispatch is carried out.

    For the \code{"formula"} method, a formula describing the structure
    of the plot, which should be of the form \code{~ x | g1 * g2 *
    \dots}, where \code{x} is a data frame or matrix.  Each of \code{g1,
    g2, \dots} must be either factors or shingles. The conditioning
    variables \code{g1, g2, \dots} may be omitted.

    For the \code{data.frame} and \code{matrix} methods, a data frame or
    matrix as appropriate.
  }
  \item{data}{
    For the \code{formula} method, an optional data frame in which
    variables in the formula (as well as \code{groups} and
    \code{subset}, if any) are to be evaluated.  By default, the
    environment where the function was called from is used.
  }

  \item{groups, subset, \dots}{ see \code{\link[lattice]{splom}}.  The
    non-standard evaluation of \code{groups} and \code{subset} only
    applies in the \code{formula} method.  Apart from arguments that
    apply to \code{splom} (many of which are only documented in
    \code{\link[lattice]{xyplot}}), additional arguments meant for
    \code{panel.hexplom} (which is an alias for
    \code{\link{panel.hexbinplot}}) may also be supplied.  Such
    arguments may include ones that control details of the hexbin
    calculations, documented in \code{\link{gplot.hexbin}}}

}
\value{

  An object of class \code{"trellis"}. The
  \code{\link[lattice:update.trellis]{update}} method can be used to
  update components of the object and the
  \code{\link[lattice:print.trellis]{print}} method (usually called by
  default) will plot it on an appropriate plotting device.

}

\seealso{
  \code{\link[lattice]{splom}}, \code{\link[lattice]{xyplot}},
  \code{\link[hexbin]{hexbinplot}},
  \code{\link[lattice]{Lattice}}, \code{\link[lattice]{panel.pairs}}
}

\author{ Deepayan Sarkar \email{Deepayan.Sarkar@R-project.org},
         Nicholas Lewin-Koh \email{nikko@hailmail.net}}

\examples{

## Simple hexplom
data(NHANES)
hexplom(~NHANES[,7:14], xbins=15)

## With colors and conditioning
hexplom(~NHANES[,9:13] | Sex, data = NHANES,
        xbins = 15, colramp = magent)

## With custom panel function
hexplom(NHANES[,9:13], xbins = 20,colramp = BTY, 
        upper.panel = panel.hexboxplot)

  
}
\keyword{hplot}

