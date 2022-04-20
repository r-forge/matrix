\name{isSymmetric-methods}
\title{Methods for Function isSymmetric in Package 'Matrix'}
\docType{methods}
%
\alias{isSymmetric-methods}
\alias{isSymmetric,symmetricMatrix-method}
\alias{isSymmetric,triangularMatrix-method}
\alias{isSymmetric,diagonalMatrix-method}
\alias{isSymmetric,unpackedMatrix-method}
\alias{isSymmetric,dgeMatrix-method}
\alias{isSymmetric,dtrMatrix-method}
\alias{isSymmetric,lgeMatrix-method}
\alias{isSymmetric,ltrMatrix-method}
\alias{isSymmetric,ngeMatrix-method}
\alias{isSymmetric,ntrMatrix-method}
\alias{isSymmetric,packedMatrix-method}
\alias{isSymmetric,dtpMatrix-method}
\alias{isSymmetric,ltpMatrix-method}
\alias{isSymmetric,ntpMatrix-method}
\alias{isSymmetric,sparseMatrix-method}
\alias{isSymmetric,dgCMatrix-method}
\alias{isSymmetric,indMatrix-method}
%
\description{
  \code{isSymmetric(M)} returns a \code{\link{logical}} indicating
  whether \code{M} is a symmetric matrix, possibly within some
  numerical tolerance.  \code{isSymmetric} is a generic function
  in \pkg{base}, which has a \link[=isSymmetric]{method} for
  traditional matrices of implicit \code{\link{class}} \code{"matrix"}.
  Methods here are defined for various proper and virtual classes
  in \pkg{Matrix}, so that \code{isSymmetric} works for all objects
  inheriting from virtual class \code{"\linkS4class{Matrix}"}.
}
\seealso{\code{\link{forceSymmetric}}, \code{\link{symmpart}}, and the
  formal class (and subclasses) \code{"\linkS4class{symmetricMatrix}"}.
}
\usage{
\S4method{isSymmetric}{symmetricMatrix}(object, \dots)
\S4method{isSymmetric}{triangularMatrix}(object, \dots)
\S4method{isSymmetric}{diagonalMatrix}(object, \dots)
\S4method{isSymmetric}{unpackedMatrix}(object, checkDN = TRUE, \dots)
\S4method{isSymmetric}{dgeMatrix}(object, tol = 100 * .Machine$double.eps, tol1 = 8 * tol, checkDN = TRUE, \dots)
\S4method{isSymmetric}{lgeMatrix}(object, checkDN = TRUE, \dots)
\S4method{isSymmetric}{ngeMatrix}(object, checkDN = TRUE, \dots)
\S4method{isSymmetric}{dtrMatrix}(object, tol = 100 * .Machine$double.eps, tol1 = 8 * tol, checkDN = TRUE, \dots)
\S4method{isSymmetric}{ltrMatrix}(object, checkDN = TRUE, \dots)
\S4method{isSymmetric}{ntrMatrix}(object, checkDN = TRUE, \dots)
\S4method{isSymmetric}{packedMatrix}(object, checkDN = TRUE, \dots)
\S4method{isSymmetric}{dtpMatrix}(object, tol = 100 * .Machine$double.eps, tol1 = 8 * tol, checkDN = TRUE, \dots)
\S4method{isSymmetric}{ltpMatrix}(object, checkDN = TRUE, \dots)
\S4method{isSymmetric}{ntpMatrix}(object, checkDN = TRUE, \dots)
\S4method{isSymmetric}{sparseMatrix}(object, tol = 100 * .Machine$double.eps, \dots)
\S4method{isSymmetric}{dgCMatrix}(object, tol = 100 * .Machine$double.eps, \dots)
\S4method{isSymmetric}{indMatrix}(object, \dots)
}
\arguments{
  \item{object}{A \code{"Matrix"}.}
  \item{tol, tol1}{Numerical tolerances allowing \emph{approximate}
    symmetry of numeric (rather than logical) matrices. See also
    \code{\link{isSymmetric.matrix}}.}
  \item{checkDN}{A logical indicating whether symmetry of the
    \code{Dimnames} slot should also be checked. \code{Dimnames}
    of the form \code{list(rn, NULL)} and \code{list(NULL, cn)}
    are considered symmetric.}
  \item{\dots}{Further arguments passed to or from methods
    (usually methods for \code{\link{all.equal}}).}
}
%TODO
%\details{
%}
\examples{
isSymmetric(Diagonal(4)) # TRUE of course
M <- Matrix(c(1,2,2,1), 2,2)
isSymmetric(M) # TRUE (*and* of formal class "dsyMatrix")
isSymmetric(as(M, "dgeMatrix")) # still symmetric, even if not "formally"
isSymmetric(triu(M)) # FALSE

## Look at implementations:
showMethods("isSymmetric", includeDefs=TRUE)# "ANY": base's S3 generic; 6 more
}
\keyword{methods}
