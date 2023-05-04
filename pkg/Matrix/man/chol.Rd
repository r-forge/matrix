\name{chol-methods}
\title{Compute the Cholesky Factor of a Matrix}
\docType{methods}
%
\keyword{algebra}
\keyword{array}
\keyword{methods}
%
\alias{chol}
\alias{chol-methods}
%
\alias{chol,diagonalMatrix-method}
\alias{chol,dgCMatrix-method}
\alias{chol,dgRMatrix-method}
\alias{chol,dgTMatrix-method}
\alias{chol,dsCMatrix-method}
\alias{chol,dsRMatrix-method}
\alias{chol,dsTMatrix-method}
\alias{chol,dspMatrix-method}
\alias{chol,dsyMatrix-method}
\alias{chol,generalMatrix-method}
\alias{chol,symmetricMatrix-method}
\alias{chol,triangularMatrix-method}
%
\description{
  Computes the upper triangular Cholesky factor of an
  \eqn{n \times n}{n-by-n} real, symmetric, positive definite
  matrix \eqn{A}, in some cases after pivoting.
  This is the factor \eqn{L'} in
  \deqn{P_{1} A P_{1}' = L L'}{P1 * A * P1' = L * L'}
  or (equivalently)
  \deqn{A = P_{1}' L L' P_{1}}{A = P1' * L * L' * P1}
  where
  \eqn{P_{1}}{P1} is a permutation matrix.

  Methods for \code{\linkS4class{denseMatrix}} are built on
  LAPACK routines \code{DPOTRF} and \code{DPPTRF}, which do
  not permute rows or columns, so that \eqn{P_{1}}{P1} is an
  identity matrix.
  
  Methods for \code{\linkS4class{sparseMatrix}} are built on
  CHOLMOD routines \code{cholmod_analyze} and \code{cholmod_factorize_p}.
}
\usage{
chol(x, \dots)
\S4method{chol}{dsyMatrix}(x, \dots)
\S4method{chol}{dspMatrix}(x, \dots)
\S4method{chol}{dsCMatrix}(x, pivot = FALSE, \dots)
\S4method{chol}{dsRMatrix}(x, pivot = FALSE, cache = TRUE, \dots)
\S4method{chol}{dsTMatrix}(x, pivot = FALSE, cache = TRUE, \dots)
}
\arguments{
  \item{x}{a \link[=is.finite]{finite}, symmetric, positive definite
    matrix or \code{\linkS4class{Matrix}} to be factorized.}
  \item{pivot}{a logical indicating if a fill-reducing
    permutation of the rows and columns of \eqn{A} should be used.
    If \code{TRUE}, then the approximate minimum degree (AMD)
    algorithm is used to select the permutation.
    If \code{FALSE}, then the order of the columns of \eqn{A}
    is preserved.}
  \item{cache}{a logical indicating if the result should be
    cached in \code{x@factors[["s[Pp]dCholesky"]]}.  Note that
    caching is experimental and that only methods for classes
    extending \code{\linkS4class{compMatrix}} will have this
    argument.}
  \item{\dots}{further arguments passed to or from methods.}
}
\value{
  A matrix or \code{\linkS4class{triangularMatrix}} representing
  the upper triangular Cholesky factor \eqn{L'}.
  The result is a traditional matrix if \code{x} is.
  For \code{x} inheriting from
  \code{\linkS4class{unpackedMatrix}},
  \code{\linkS4class{packedMatrix}}, and
  \code{\linkS4class{sparseMatrix}},
  the result is an object of class
  \code{\linkS4class{Cholesky}},
  \code{\linkS4class{pCholesky}}, and
  \code{\linkS4class{dtCMatrix}},
  respectively.
}
\details{
  For sparse \code{x}, \code{chol(x, pivot = value)} uses
  \code{\link{Cholesky}(x, perm = value, LDL = FALSE)} under the hood.
  If you must know the permutation \eqn{P_{1}}{P1} in addition to
  the Cholesky factor \eqn{L} (or in this case the transpose \eqn{L'}),
  then call \code{\link{Cholesky}} directly, as the result of
  \code{chol(x, pivot = TRUE)} specifies \eqn{L'} but not \eqn{P_{1}}{P1}.
}
\seealso{
  The default method from \pkg{base}, \code{\link[base]{chol}},
  called for traditional matrices \code{x}.

  Generic function \code{\link{Cholesky}}, for more flexibility
  notably when computing the Cholesky \emph{factorization} and
  not only the \emph{factor} \eqn{L'}.
}
\references{
  The LAPACK source code, including documentation; see
  \url{https://netlib.org/lapack/double/dpotrf.f} and
  \url{https://netlib.org/lapack/double/dpptrf.f}.

  Davis, T. A. (2006).
  \emph{Direct Methods for Sparse Linear Systems},
  Society for Industrial and Applied Mathematics.
  
  Amestoy, P. R., Davis, T. A., & Duff, I. S. (2004).
  Algorithm 837: AMD, an approximate minimum degree ordering algorithm.
  \emph{ACM Transactions on Mathematical Software},
  \emph{17}(4), 886-905.
  \doi{10.1145/1024074.1024081}
}
\examples{
showMethods("chol", inherited = FALSE)

sy2 <- new("dsyMatrix", Dim = as.integer(c(2,2)), x = c(14, NA,32,77))
(c2 <- chol(sy2))#-> "Cholesky" matrix
stopifnot(all.equal(c2, chol(as(sy2, "dpoMatrix")), tolerance= 1e-13))
str(c2)

## An example where chol() can't work
(sy3 <- new("dsyMatrix", Dim = as.integer(c(2,2)), x = c(14, -1, 2, -7)))
try(chol(sy3)) # error, since it is not positive definite

## A sparse example --- exemplifying 'pivot'
(mm <- toeplitz(as(c(10, 0, 1, 0, 3), "sparseVector"))) # 5 x 5
(R <- chol(mm)) ## default:  pivot = FALSE
R2 <- chol(mm, pivot=FALSE)
stopifnot( identical(R, R2), all.equal(crossprod(R), mm) )
(R. <- chol(mm, pivot=TRUE))# nice band structure,
## but of course crossprod(R.) is *NOT* equal to mm
## --> see Cholesky() and its examples, for the pivot structure & factorization
stopifnot(all.equal(sqrt(det(mm)), det(R)),
          all.equal(prod(diag(R)), det(R)),
          all.equal(prod(diag(R.)), det(R)))

## a second, even sparser example:
(M2 <- toeplitz(as(c(1,.5, rep(0,12), -.1), "sparseVector")))
c2 <- chol(M2)
C2 <- chol(M2, pivot=TRUE)
## For the experts, check the caching of the factorizations:
ff <- M2@factors[["spdCholesky"]]
FF <- M2@factors[["sPdCholesky"]]
L1 <- as(ff, "Matrix")# pivot=FALSE: no perm.
L2 <- as(FF, "Matrix"); P2 <- as(FF, "pMatrix")
stopifnot(identical(t(L1), c2),
          all.equal(t(L2), C2, tolerance=0),#-- why not identical()?
          all.equal(M2, tcrossprod(L1)),             # M = LL'
          all.equal(M2, crossprod(crossprod(L2, P2)))# M = P'L L'P
         )
}
