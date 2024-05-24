\name{forceSymmetric-methods}
\title{Symmetrize a Square Matrix}
%
\docType{methods}
\keyword{array}
\keyword{methods}
%
\alias{forceSymmetric}
\alias{forceSymmetric-methods}
%
\alias{forceSymmetric,CsparseMatrix-method}
\alias{forceSymmetric,RsparseMatrix-method}
\alias{forceSymmetric,TsparseMatrix-method}
\alias{forceSymmetric,denseMatrix-method}
\alias{forceSymmetric,diagonalMatrix-method}
\alias{forceSymmetric,indMatrix-method}
\alias{forceSymmetric,matrix-method}
%
\description{
Returns the symmetric matrix that results from reflecting one triangle
of a square matrix onto the other.
}
\usage{
forceSymmetric(x, \dots)

\S4method{forceSymmetric}{denseMatrix}(x, uplo = NULL, trans = "C", \dots)
\S4method{forceSymmetric}{CsparseMatrix}(x, uplo = NULL, trans = "C", \dots)
\S4method{forceSymmetric}{RsparseMatrix}(x, uplo = NULL, trans = "C", \dots)
\S4method{forceSymmetric}{TsparseMatrix}(x, uplo = NULL, trans = "C", \dots)
\S4method{forceSymmetric}{diagonalMatrix}(x, uplo = "U", trans = "C", \dots)
\S4method{forceSymmetric}{indMatrix}(x, \dots)
\S4method{forceSymmetric}{matrix}(x, uplo = "U", trans = "C", \dots)
}
\arguments{
\item{x}{a square matrix or \code{\linkS4class{Matrix}}.}
\item{uplo}{a character string, \code{"U"} or \code{"L"}, indicating
  which of the upper and lower triangles should be reflected.
  Where it is the default, \code{NULL} is accepted and equivalent
  to \code{"U"}, unless \code{x} has an \code{uplo} slot, in which
  case that is used instead.}
\item{trans}{a character string, \code{"C"} or \code{"T"},
  indicating either the conjugate transpose or transpose,
  i.e., whether the result is formally Hermitian or symmetric.
  It matters only in the complex case.  The default is \code{"C"},
  for consistency with methods for \code{\link{isSymmetric}}.}
\item{\dots}{optional arguments passed from or to other methods.}
}
\value{
An object inheriting from virtual class
\code{\linkS4class{symmetricMatrix}},
representing a Hermitian or symmetric matrix.
}
\details{
\code{forceSymmetric(x)} can be understood as a coercion to virtual
class \code{\linkS4class{symmetricMatrix}}.
However, unlike \code{as(x, "symmetricMatrix")},
methods for \code{forceSymmetric} do \emph{not} require that \code{x}
is already approximately symmetric.
}
\seealso{
\code{\link{isSymmetric}}, \code{\link{symmpart}}, \code{\link{skewpart}}.
}
\examples{
 ## Hilbert matrix
 i <- 1:6
 h6 <- 1/outer(i - 1L, i, "+")
 sd <- sqrt(diag(h6))
 hh <- t(h6/sd)/sd # theoretically symmetric
 isSymmetric(hh, tol=0) # FALSE; hence
 try( as(hh, "symmetricMatrix") ) # fails, but this works fine:
 H6 <- forceSymmetric(hh)

 ## result can be pretty surprising:
 (M <- Matrix(1:36, 6))
 forceSymmetric(M) # symmetric, hence very different in lower triangle
 (tm <- tril(M))
 forceSymmetric(tm)
}
