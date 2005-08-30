\name{simulate}
\docType{genericFunction}
\alias{simulate}
\alias{simulate,lmer-method}
\title{Simulate responses from a fitted model}
\description{
  This generic function generates a sample of responses from a fitted
  model object using the parameter estimates from the fit.  Typically
  these simulated responses are used in a parametric bootstrap sample of
  a test statistic.
}
\usage{
simulate(object, n, \dots)
}
\arguments{
  \item{object}{An object of a suitable class - usually an
    \code{\link{lmer-class}{lmer}} object. 
  }
  \item{n}{integer - number of samples to simulate. Defaults to 1.}
  \item{\dots}{Some methods for this generic function may take
    additional, optional arguments.}
}
\value{
  A matrix of dimension \code{(nobs, n)} where \code{nobs} is the number
  of observations in the fitted model \code{object}.
}
\section{Methods}{
  \describe{
    \item{object = "lmer"}{simulate responses from the fitted model \code{"object"}.}
  }
}
\examples{
(fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0 + Days|Subject), sleepstudy))
sim2 <- simulate(fm2, n = 4)
head(sim2)
tail(sim2)
}
\keyword{methods}
\keyword{datagen}
