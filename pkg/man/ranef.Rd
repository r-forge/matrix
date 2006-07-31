% $Id: random.effects.Rd,v 1.1 2002/11/27 22:49:45 bates Exp $
\name{ranef}
\docType{genericFunction}
\docType{methods}
\alias{ranef}
\alias{random.effects}
\alias{ranef-methods}
\alias{ranef,ANY-method}
\alias{ranef,lmer-method}
\alias{ranef,mer-method}
\title{Extract Random Effects}
\usage{
ranef(object, \dots)
\S4method{ranef}{lmer}(object, postVar, \dots)
random.effects(object, \dots)   # deprecated
}
\description{
  A generic function to extract, and optionally accumulate, the random effects.
}
\arguments{
 \item{object}{an object of a class from which random effects
   estimates can be extracted.}
 \item{postVar}{an optional logical argument indicating if the
   conditional variance covariance matrices, also called the "posterior
   variances", of the random effects should be included.  Default is
   \code{FALSE}.}
 \item{\dots}{some methods for this generic function require additional
   arguments.}
}
\value{
  A list of data frames, one for each grouping factor for the random
  effects.  The number of rows in the data frame is the number of levels
  of the grouping factor.  The number of columns is the dimension of the
  random effect associated with each level of the factor.

  If \code{postVar} is \code{TRUE} each of the data frames has an
  attribute called \code{"postVar"} which is a three-dimensional array
  with symmetric faces.
}
\examples{
fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy)
ranef(fm2)
str(ranef(fm2, postVar = TRUE))
}
\keyword{models}
\keyword{methods}
