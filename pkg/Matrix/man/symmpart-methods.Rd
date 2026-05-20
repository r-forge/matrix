\name{symmpart-methods}
%%
\alias{symmpart}
\alias{symmpart-methods}
\alias{symmpart,CsparseMatrix-method}
\alias{symmpart,RsparseMatrix-method}
\alias{symmpart,TsparseMatrix-method}
\alias{symmpart,denseMatrix-method}
\alias{symmpart,diagonalMatrix-method}
\alias{symmpart,indMatrix-method}
\alias{symmpart,matrix-method}
%%
\alias{skewpart}
\alias{skewpart-methods}
\alias{skewpart,CsparseMatrix-method}
\alias{skewpart,RsparseMatrix-method}
\alias{skewpart,TsparseMatrix-method}
\alias{skewpart,denseMatrix-method}
\alias{skewpart,diagonalMatrix-method}
\alias{skewpart,indMatrix-method}
\alias{skewpart,matrix-method}
%%
\title{Symmetric and Skew-Symmetric Parts of a Square Matrix}
\description{
  Compute the symmetric (Hermitian) and skew-symmetric (skew-Hermitian)
  parts of a square matrix \code{x}, defined as \code{(x + op(x))/2} and
  \code{(x - op(x))/2}, respectively, where \code{op} denotes matrix
  transpose (conjugate transpose).
}
\usage{
symmpart(x, \dots)
skewpart(x, \dots)

\S4method{symmpart}{denseMatrix}(x, uplo = "U", trans = "C", \dots)
\S4method{symmpart}{CsparseMatrix}(x, uplo = "U", trans = "C", \dots)
\S4method{symmpart}{RsparseMatrix}(x, uplo = "U", trans = "C", \dots)
\S4method{symmpart}{TsparseMatrix}(x, uplo = "U", trans = "C", \dots)
\S4method{symmpart}{diagonalMatrix}(x,             trans = "C", \dots)
\S4method{symmpart}{indMatrix}(x,                          \dots)
\S4method{symmpart}{matrix}(x,             trans = "C", \dots)

\S4method{skewpart}{denseMatrix}(x, trans = "C", \dots)
\S4method{skewpart}{CsparseMatrix}(x, trans = "C", \dots)
\S4method{skewpart}{RsparseMatrix}(x, trans = "C", \dots)
\S4method{skewpart}{TsparseMatrix}(x, trans = "C", \dots)
\S4method{skewpart}{diagonalMatrix}(x, trans = "C", \dots)
\S4method{skewpart}{indMatrix}(x,              \dots)
\S4method{skewpart}{matrix}(x, trans = "C", \dots)
}
\arguments{
  \item{x}{an \R{} object representing a square matrix.  Package
    \pkg{Matrix} provides methods for \code{\linkS4class{Matrix}}
    objects and traditional matrices.}
  \item{uplo}{a character string, \code{"U"} or \code{"L"}, indicating
    which of the upper and lower triangles of the symmetric (Hermitian)
    part is stored in the result.}
  \item{trans}{a character string, \code{"C"} or \code{"T"}, indicating
    either conjugate transpose or transpose, i.e., whether the result is
    (skew-)Hermitian or (skew-)symmetric.  It matters only for \code{x}
    with complex data.  The default value is \code{"C"} for consistency
    with methods for \code{\link{isSymmetric}}.}
  \item{\dots}{optional arguments passed from or to other methods.}
}
\value{
  An \R{} object representing a symmetric (Hermitian) or skew-symmetric
  (skew-Hermitian) matrix, inheriting from virtual class
  \code{\linkS4class{Matrix}} if \code{x} does.  The dimension names are
  the result of symmetrizing \code{dimnames(x)} giving precedence to the
  column names unless those are \code{NULL} and the row names are
  non-\code{NULL}.
}
\seealso{
  \code{\linkS4methods{isSymmetric}},
  \code{\linkS4methods{forceSymmetric}}.
}
\examples{
(X <- new("zgeMatrix",
          Dim = c(3L, 3L),
          Dimnames = list(A = paste0("a", 1:3), B = paste0("b", 1:3)),
          x = 0:8 * rep(c(1, -1i), length.out = 9L)))
(Xc0 <- symmpart(X))              #      Hermitian
(Xc1 <- skewpart(X))              # skew-Hermitian
(Xt0 <- symmpart(X, trans = "T")) #      symmetric
(Xt1 <- skewpart(X, trans = "T")) # skew-symmetric
stopifnot(isSymmetric(Xc0),
          isSymmetric(Xc1 * 1i),
          isSymmetric(Xt0, trans = "T"),
          isSymmetric(Xt1 * ifelse(upper.tri(Xt1), 1, -1), trans = "T"),
          all(Xc0 + Xc1 == X),
          all(Xt0 + Xt1 == X),
          all(Im(diag(Xc0)) == 0),
          all(Re(diag(Xc1)) == 0),
          all(diag(Xt1) == 0i),
          max(abs(Im(eigen(Xc0)$values))) < 16 * .Machine$double.eps,
          max(abs(Re(eigen(Xc1)$values))) < 16 * .Machine$double.eps,
          identical(lapply(list(Xc0, Xc1, Xt0, Xt1), dimnames),
                    rep(list(dimnames(X)[c(2L, 2L)]), 4L)))
}
\docType{methods}
\keyword{algebra}
\keyword{arith}
\keyword{array}
\keyword{methods}
