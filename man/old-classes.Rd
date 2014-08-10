\name{old-classes}
\title{Class "unit" and "viewport" as S4 classes}
%
\docType{class}
\alias{unit-class}
\alias{viewport-class}
%
\description{Package "hexbin" now uses S4 classes throughout and hence
  needs to \code{\link[methods]{setOldClass}} both \code{"unit"} and
  \code{"viewport"} (which are S3 classes from the \pkg{grid} package),
  in order to be able to use those in slots of its own classes.
}
\section{Objects from the Class}{A virtual Class: No objects may be
  created from it.}
\section{Extends}{
  Class \code{"oldClass"}, directly.
}
\section{Methods}{
  No methods defined with class "unit" in the signature.
}
% \seealso{
%               add link to grid ??
% }
\keyword{classes}
