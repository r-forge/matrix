\name{forceCanonical-methods}
\title{Canonicalize an Object}
%
\docType{methods}
\keyword{array}
\keyword{programming}
\keyword{methods}
%
\alias{forceCanonical}
\alias{forceCanonical-methods}
%
\alias{forceCanonical,CsparseMatrix-method}
\alias{forceCanonical,RsparseMatrix-method}
\alias{forceCanonical,TsparseMatrix-method}
\alias{forceCanonical,denseMatrix-method}
\alias{forceCanonical,diagonalMatrix-method}
\alias{forceCanonical,indMatrix-method}
\alias{forceCanonical,sparseVector-method}
%
\description{
  Obtain a \dQuote{canonical} representation of an \R{} object as
  defined by corresponding methods for \code{\link{isCanonical}}.
}
\usage{
forceCanonical(x, \dots)
}
\arguments{
  \item{x}{an \R{} object.}
  \item{\dots}{optional arguments passed to methods.}
}
\details{
  Method writers should ensure that canonicalization is idempotent.
  In other words, \code{f(f(x))} and \code{f(x)} should be identical
  if \code{f = forceCanonical}.
}
\value{
  A canonical object, typically of class \code{class(x)} or a
  superclass, such that \code{isCanonical(.)} is \code{TRUE}.
}
\seealso{\code{\link{isCanonical}}.}
\examples{
utils::example(isCanonical, package = "Matrix")
stopifnot(identical(unname(obj), lapply(unname(obj.), forceCanonical)))
}
