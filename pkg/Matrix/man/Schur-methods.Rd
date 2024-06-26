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
\alias{Schur,denseMatrix-method}
\alias{Schur,matrix-method}
\alias{Schur,sparseMatrix-method}
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
  An object representing the factorization, inheriting
  from virtual class \code{SchurFactorization}.
  Currently, the specific class is always \code{Schur}.
}
\seealso{
  Class \code{Schur} and its methods.

  Class \code{\linkS4class{dgeMatrix}}.

  Generic functions \code{\link{expand1}} and \code{\link{expand2}},
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
\dontshow{ % for R_DEFAULT_PACKAGES=NULL
library(stats, pos = "package:base", verbose = FALSE)
library(utils, pos = "package:base", verbose = FALSE)
}
showMethods("Schur", inherited = FALSE)
set.seed(0)

Schur(Hilbert(9L)) # real eigenvalues

n <- 5L
(A <- Matrix(round(rnorm(n * n, sd = 100)), n, n))
(sch.A <- Schur(A)) # complex eigenvalues

## A ~ Q T Q' in floating point
str(e.sch.A <- expand2(sch.A), max.level = 2L)
stopifnot(all.equal(A, Reduce(`\%*\%`, e.sch.A)))

(e1 <- eigen(matrix(sch.A@x, n, n), only.values = TRUE)$values)
(e2 <- eigen(           A         , only.values = TRUE)$values)
(e3 <- sch.A@values)

stopifnot(exprs = {
    all.equal(e1, e2, tolerance = 1e-13)
    all.equal(e1, e3[order(Mod(e3), decreasing = TRUE)], tolerance = 1e-13)
    identical(Schur(A, vectors = FALSE),
              `slot<-`(sch.A, "vectors", TRUE, double(0L)))
})
}
