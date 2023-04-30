\name{BunchKaufman-methods}
\title{Methods for Bunch-Kaufman Factorization}
%
\docType{methods}
\keyword{algebra}
\keyword{array}
\keyword{methods}
%
\alias{BunchKaufman}
\alias{BunchKaufman-methods}
%
\alias{BunchKaufman,dspMatrix-method}
\alias{BunchKaufman,dsyMatrix-method}
\alias{BunchKaufman,matrix-method}
%
\description{
  The Bunch-Kaufman factorization of a symmetric matrix \eqn{A}
  is a factorization of the form \eqn{A = U' D U} or \eqn{A = L D L'}
  where \eqn{U} and \eqn{L} are themselves factored as products of
  permutation matrices and unit upper or lower triangular matrices
  and where \eqn{D} is a block diagonal matrix composed of 1-by-1 and
  2-by-2 diagonal blocks.
  
  This factorization generalizes the pivoted Cholesky factorization
  to symmetric matrices which may not be positive semidefinite.
}
\usage{
\S4method{BunchKaufman}{dsyMatrix}(x, \dots)
\S4method{BunchKaufman}{dspMatrix}(x, \dots)
\S4method{BunchKaufman}{matrix}(x, uplo = "U", \dots)
}
\arguments{
  \item{x}{a square matrix.}
  \item{uplo}{a string, either \code{"U"} or \code{"L"}, indicating
    which triangle (upper or lower) of \code{x} should be used to
    compute the factorization.}
  \item{\dots}{further optional arguments passed to methods.}
}
\section{Methods}{
  \describe{
    \item{\code{x = "dsyMatrix"}}{using LAPACK routine \code{DSYTRF}.}
    \item{\code{x = "dspMatrix"}}{using LAPACK routine \code{DSPTRF}.}
    \item{\code{x = "matrix"}}{using LAPACK routine \code{DSYTRF},
      where \code{x} must have type \link{double} and must be square
      but need not be symmetric.}
  }
}
\value{
  An object inheriting from class \code{\linkS4class{BunchKaufman}}
  or \code{\linkS4class{pBunchKaufman}}, the latter if and only if
  \code{x} inherits from \code{\linkS4class{dspMatrix}}.
}
\references{
  The original LAPACK source code, including documentation; see
  \url{https://netlib.org/lapack/double/dsytrf.f} and
  \url{https://netlib.org/lapack/double/dsptrf.f}.
}
\seealso{
  Classes \code{\linkS4class{dsyMatrix}} and
  \code{\linkS4class{BunchKaufman}}.
  
  Related are the LU factorization defined for general matrices
  and computed by \code{\link{lu}}, and the Cholesky factorization
  defined for symmetric positive definite matrices and computed by
  \code{\link{chol}} or (in the sparse case) \code{\link{Cholesky}}.
}
\examples{
data(CAex)
dim(CAex)
isSymmetric(CAex) # TRUE

S <- as(CAex, "symmetricMatrix") # dsCMatrix

## No methods for sparseMatrix yet:
if(FALSE)
   BunchKaufman(S)
## But two for *dense*Matrix:
 bk <- BunchKaufman(as(S, "unpackedMatrix"))
pbk <- BunchKaufman(as(S,   "packedMatrix"))
pbk

## See the underlying LAPACK representation
 bk.M <- as( bk, "dtrMatrix")
pbk.M <- as(pbk, "dtpMatrix")
pbk.M[1:6, 1:6]

## See how sparse it is, too
bk.M.s <- as(bk.M, "sparseMatrix")
image(bk.M.s)

## An upper triangular 72-by-72 matrix with only 144 nonzero entries
str(bk.M.s)

stopifnot(is(bk.M,   "triangularMatrix"),
          is(bk.M,     "unpackedMatrix"),
          is(bk.M.s, "triangularMatrix"),
          is(bk.M.s,    "CsparseMatrix"))
}
