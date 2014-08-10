\name{ColorRamps}
\title{Color Ramps on Perceptually Linear Scales}
\alias{ColorRamps}
\alias{LinGray}
\alias{BTC}
\alias{BTY}
\alias{LinOCS}
\alias{heat.ob}
\alias{magent}
\alias{plinrain}
\description{
  Functions for returning colors on perceptually linear scales,
  where steps correspond to \sQuote{just detectable differences}.
}
\usage{
LinGray (n, beg=1, end=92)
BTC     (n, beg=1, end=256)
LinOCS  (n, beg=1, end=256)
heat.ob (n, beg=1, end=256)
magent  (n, beg=1, end=256)
plinrain(n, beg=1, end=256)
}
\arguments{
  \item{n}{number of colors to return from the ramp}
  \item{beg}{begining of ramp, integer from 1-255}
  \item{end}{end of ramp, integer from 1-255}
}
\value{
  returns an array of colors
}
\details{
  Several precalulated color ramps, that are on a perceptually linear
  color scale.  A perceptually linear color scale is a scale where each
  jump corresponds to a \dQuote{just detectable difference} in color and the
  scale is percieved as linear by the human eye (emprically determined).

  When using the ramps, if \code{beg} is less than \code{end} the ramp
  will be reversed.
}
\references{
  Haim Levkowitz (1997)
  \emph{Color Theory and Modeling for Computer Graphics,
    Visualization, and Multimedia Applications}.
  Kluwer Academic Publishers, Boston/London/Dordrecht.
  \url{http://www.cs.uml.edu/~haim/ColorCenter/}
}
\seealso{
  \code{\link[grDevices:palettes]{rainbow}},
  \code{\link[grDevices:palettes]{terrain.colors}},
  \code{\link[grDevices]{rgb}},  
  \code{\link[grDevices]{hsv}}
}
\examples{
h <- hexbin(rnorm(10000),rnorm(10000))
plot(h, colramp= BTY)
## looks better if you shave the tails:
plot(h, colramp= function(n){LinOCS(n,beg=15,end=225)})
}
\author{Nicholas Lewin-Koh}
\keyword{color}
