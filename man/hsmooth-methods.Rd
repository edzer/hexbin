\name{hsmooth-methods}
\docType{methods}% + generic -- still use this doctype ?
\alias{hsmooth}% generic
\alias{hsmooth-methods}
\alias{hsmooth,hexbin-method}
\title{Hexagon Bin Smoothing: Generic hsmooth() and Methods}
\description{
  Methods for the generic function \code{hsmooth} in package
  \pkg{hexbin}:
  There is currently only the one for \code{\link{hexbin}} objects.
}
\usage{
\S4method{hsmooth}{hexbin}(bin, wts)
}
\arguments{
  \item{bin}{a \code{\link{hexbin}} object, or an extension such as
    \code{\link{erodebin-class}}.}
  \item{wts}{weights vector, see \code{\link{smooth.hexbin}}}
}
\section{Methods}{
  \describe{
    \item{bin = "hexbin"}{is just the \code{\link{smooth.hexbin}}
      function (for back compatibility); see its documentation, also for
      examples.}
  }
}
\keyword{methods}

