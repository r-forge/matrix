\name{Schur-methods}
\title{Methods for Schur Factorization}
%
\docType{methods}
\keyword{algebra}
\keyword{array}
\keyword{methods}
%
\alias{Schur}
\alias{Schur-methods}
%
\alias{Schur,Matrix,missing-method}
\alias{Schur,matrix,missing-method}
\alias{Schur,dgeMatrix,logical-method}
\alias{Schur,diagonalMatrix,logical-method}
\alias{Schur,dsyMatrix,logical-method}
\alias{Schur,generalMatrix,logical-method}
\alias{Schur,matrix,logical-method}
\alias{Schur,symmetricMatrix,logical-method}
\alias{Schur,triangularMatrix,logical-method}
%
\description{
  Computes the Schur factorization of an \eqn{n \times n}{n-by-n}
  real matrix \eqn{A}, which has the general form 
  \deqn{A = Q T Q'}{A = Q * T * Q'}
  where
  \eqn{Q} is an orthogonal matrix and
  \eqn{T} is a block upper triangular matrix with
  \eqn{1 \times 1}{1-by-1} and \eqn{2 \times 2}{2-by-2} diagonal blocks
  specifying the real and complex conjugate eigenvalues of \eqn{A}.
  The column vectors of \eqn{Q} are the Schur vectors of \eqn{A},
  and \eqn{T} is the Schur form of \eqn{A}.
  
  Methods are built on LAPACK routine \code{dgees}.
}
\usage{
Schur(x, vectors = TRUE, \dots)
}
\arguments{
  \item{x}{a \link[=is.finite]{finite} square matrix or
    \code{\linkS4class{Matrix}} to be factorized.}
  \item{vectors}{a logical.  If \code{TRUE} (the default),
    then Schur vectors are computed in addition to the Schur form.}
  \item{\dots}{further arguments passed to or from methods.}
}
\value{
  If \code{vectors = TRUE} and \code{x} is a
  \code{\linkS4class{Matrix}}, then the result is an object
  of class \code{\linkS4class{Schur}}.

  If \code{vectors = TRUE} and \code{x} is a tradtional matrix,
  then the result is a named list containing the \code{Q}, \code{T},
  and \code{EValues} slots of the \code{\linkS4class{Schur}} object.

  If \code{vectors = FALSE}, then the result is the same named list
  but without \code{Q}.
}
\seealso{
  Class \code{\linkS4class{Schur}} and its methods.

  Class \code{\linkS4class{dgeMatrix}}.

  Generic function \code{\link{expand2}},
  for constructing matrix factors from the result.
  
  Generic functions \code{\link{Cholesky}}, \code{\link{BunchKaufman}},
  \code{\link{lu}}, and \code{\link{qr}},
  for computing other factorizations.
}
\references{
  The LAPACK source code, including documentation; see
  \url{https://netlib.org/lapack/double/dgees.f}.

  Golub, G. H., & Van Loan, C. F. (2013).
  \emph{Matrix computations} (4th ed.).
  Johns Hopkins University Press.
  \doi{10.56021/9781421407944}
}
\examples{
showMethods("Schur", inherited = FALSE)
set.seed(0)

Schur(Hilbert(9L)) # real eigenvalues

A <- Matrix(round(rnorm(25L, sd = 100)), 5L, 5L)
(Sch.A <- Schur(A)) # complex eigenvalues

## A ~ Q T Q' in floating point
e.Sch.A <- expand2(Sch.A)
names(e.Sch.A)
stopifnot(all.equal(A, Reduce(`\%*\%`, e.Sch.A)))

(e1 <- eigen(Sch.A@T, only.values = TRUE)$values)
(e2 <- eigen(    A  , only.values = TRUE)$values)
(e3 <- Sch.A@EValues)

stopifnot(exprs = {
    all.equal(e1, e2, tolerance = 1e-13)
    all.equal(e1, e3[order(Mod(e3), decreasing = TRUE)], tolerance = 1e-13) 
    identical(Schur(A, vectors = FALSE),
              list(T = Sch.A@T, EValues = e3))    
    identical(Schur(as(A, "matrix")),
              list(Q = as(Sch.A@Q, "matrix"),
                   T = as(Sch.A@T, "matrix"), EValues = e3))
})
}
