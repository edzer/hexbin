\name{NHANES}
\alias{NHANES}
\docType{data}
\title{NHANES Data : National Health and Nutrition Examination Survey}
\usage{data(NHANES)}
\description{
  This is a somewhat large interesting dataset, a data frame of 15
  variables (columns) on 9575 persons (rows).
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{Cancer.Incidence}{binary factor with levels \code{No} and \code{Yes}.}
    \item{Cancer.Death}{binary factor with levels \code{No} and \code{Yes}.}
    \item{Age}{numeric vector giving age of the person in years.}
    \item{Smoke}{a factor with levels \code{Current}, \code{Past},
      \code{Nonsmoker}, and \code{Unknown}.}
    \item{Ed}{numeric vector of \eqn{\{0,1\}} codes giving the education level.}
    \item{Race}{numeric vector of \eqn{\{0,1\}} codes giving the
      person's race.%% FIXME : 0 = ?  1 = ?
    }
    \item{Weight}{numeric vector giving the weight in kilograms}
    \item{BMI}{numeric vector giving Body Mass Index, i.e.,
      \code{Weight/Height^2} where Height is in meters, and missings
      (61\% !) are coded as \code{0} originally.}%% rather FIXME?
    \item{Diet.Iron}{numeric giving Dietary iron.}
    \item{Albumin}{numeric giving albumin level in g/l.}
    \item{Serum.Iron}{numeric giving Serum iron in \eqn{\mu}{u}g/l.}
    \item{TIBC}{numeric giving Total Iron Binding Capacity in \eqn{\mu}{u}g/l.}
    \item{Transferin}{numeric giving Transferin Saturation which is just
      \code{100*serum.iron/TIBC}.}
    \item{Hemoglobin}{numeric giving Hemoglobin level.}
    \item{Sex}{a factor with levels \code{F} (female) and \code{M} (male).}
  }
}
\examples{
data(NHANES)
summary(NHANES)
## Missing Data overview :
nNA <- sapply(NHANES, function(x)sum(is.na(x)))
cbind(nNA[nNA > 0])
# Which are just these 6 :
\dontrun{
Diet.Iron         141
Albumin           252
Serum.Iron       1008
TIBC              853
Transferin       1019
Hemoglobin        759
}%dont
}
\keyword{datasets}
