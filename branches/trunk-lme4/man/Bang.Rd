\name{Bang}
\alias{Bang}
\docType{data}
\title{Bangladesh fertility survey}
\description{
  These data come from the 1988 Bangladesh Fertility Survey.
}
\usage{data(Bang)}
\format{
  A data frame with 1934 observations on the following 6 variables.
  \describe{
    \item{woman}{Identifying code for each woman - a factor}
    \item{district}{Identifying code for each district - a factor}
    \item{use}{Contraceptive use at time of survey}
    \item{livch}{Number of living children at time of survey - an
      ordered factor.  Levels are \code{0}, \code{1}, \code{2}, \code{3+}}
    \item{age}{Age of woman at time of survey (in years), centred around
      mean.}
    \item{urban}{Type of region of residence - a factor.  Levels are
      \code{urban} and \code{rural}}
  }
}
%\details{}
\source{
  \url{http://multilevel.ioe.ac.uk/softrev/bang.html}
}
\references{
  Steele, F., Diamond, I. And Amin, S. (1996). Immunization uptake in
  rural Bangladesh: a multilevel analysis. \emph{Journal of the Royal
    Statistical Society, Series A} (159): 289-299.
}
\examples{
data(Bang)
str(Bang)
summary(Bang)
library(lme4)
system.time(fm1 <- GLMM(use ~ urban+age+livch,
   family = binomial, data = Bang, random = ~ 1|district))
summary(fm1)
system.time(fm2 <- GLMM(use ~ urban+age+livch,
   family = binomial, data = Bang, random = ~ 1|district, method = "Laplace"))
summary(fm2)
\dontrun{
system.time(fm2 <- GLMM(use ~ urban+age+livch,
   family = binomial, data = Bang, random = ~ urban|district,
   niter = 50, verbose = TRUE))
summary(fm2)
}
}
\keyword{datasets}
