\name{symmpart-methods}
\title{Symmetric and Skew-Symmetric Parts of a Square Matrix}
%
\docType{methods}
\keyword{algebra}
\keyword{arith}
\keyword{array}
\keyword{methods}
%
\alias{symmpart}
\alias{symmpart-methods}
\alias{skewpart}
\alias{skewpart-methods}
%
\alias{symmpart,CsparseMatrix-method}
\alias{symmpart,RsparseMatrix-method}
\alias{symmpart,TsparseMatrix-method}
\alias{symmpart,denseMatrix-method}
\alias{symmpart,diagonalMatrix-method}
\alias{symmpart,indMatrix-method}
\alias{symmpart,matrix-method}
%
\alias{skewpart,CsparseMatrix-method}
\alias{skewpart,RsparseMatrix-method}
\alias{skewpart,TsparseMatrix-method}
\alias{skewpart,denseMatrix-method}
\alias{skewpart,diagonalMatrix-method}
\alias{skewpart,indMatrix-method}
\alias{skewpart,matrix-method}
%
\description{
\code{symmpart} and \code{skewpart} compute the symmetric and skew
symmetric parts \code{(x + t(x))/2} and \code{(x - t(x))/2} of a
square matrix \code{x}.  \code{symmpart(x) + skewpart(x)} is equal
to \code{x} for all finite square matrices \code{x}.
}
\usage{
symmpart(x, \dots)
skewpart(x, \dots)

\S4method{symmpart}{denseMatrix}(x, uplo = "U", trans = "C", \dots)
\S4method{symmpart}{CsparseMatrix}(x, uplo = "U", trans = "C", \dots)
\S4method{symmpart}{RsparseMatrix}(x, uplo = "U", trans = "C", \dots)
\S4method{symmpart}{TsparseMatrix}(x, uplo = "U", trans = "C", \dots)
\S4method{symmpart}{diagonalMatrix}(x, trans = "C", \dots)
\S4method{symmpart}{indMatrix}(x, \dots)
\S4method{symmpart}{matrix}(x, trans = "C", \dots)

## Methods for 'skewpart' are parallel, but without 'uplo'
}
\arguments{
\item{x}{a square matrix or \code{\linkS4class{Matrix}}.}
\item{uplo}{a character string, \code{"U"} or \code{"L"}, indicating
  which of the upper and lower triangles of the symmetric part is
  stored in the result.}
\item{trans}{a character string, \code{"C"} or \code{"T"},
  indicating either the conjugate transpose or transpose,
  i.e., whether the result is (skew) Hermitian or symmetric.
  It matters only in the complex case.  The default is \code{"C"},
  for consistency with methods for \code{\link{isSymmetric}}.}
\item{\dots}{optional arguments passed from or to other methods.}
}
\value{
\code{symmpart(x)} returns a symmetric or Hermitian matrix,
inheriting from \code{\linkS4class{symmetricMatrix}}
or \code{\linkS4class{diagonalMatrix}} if \code{x}
inherits from \code{Matrix}.
  
\code{skewpart(x)} returns a skew-symmetric or skew-Hermitian matrix,
inheriting from \code{\linkS4class{generalMatrix}},
\code{\linkS4class{symmetricMatrix}} or
\code{\linkS4class{diagonalMatrix}} if \code{x}
inherits from \code{Matrix}.
}
\seealso{
\code{\link{isSymmetric}}, \code{\link{forceSymmetric}}.
}
\examples{
m <- Matrix(1:4, 2,2)
symmpart(m)
skewpart(m)

stopifnot(all(m == symmpart(m) + skewpart(m)))

dn <- dimnames(m) <- list(row = c("r1", "r2"), col = c("var.1", "var.2"))
stopifnot(all(m == symmpart(m) + skewpart(m)))
colnames(m) <- NULL
stopifnot(all(m == symmpart(m) + skewpart(m)))
dimnames(m) <- unname(dn)
stopifnot(all(m == symmpart(m) + skewpart(m)))


## investigate the current methods:
showMethods(skewpart, include = TRUE)
}
