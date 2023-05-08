\name{lu-methods}
\title{Methods for LU Factorization}
%
\docType{methods}
\keyword{algebra}
\keyword{array}
\keyword{methods}
%
\alias{lu}
\alias{lu-methods}
%
\alias{lu,denseMatrix-method}
\alias{lu,diagonalMatrix-method}
\alias{lu,dgCMatrix-method}
\alias{lu,dgRMatrix-method}
\alias{lu,dgTMatrix-method}
\alias{lu,dgeMatrix-method}
\alias{lu,dsCMatrix-method}
\alias{lu,dsRMatrix-method}
\alias{lu,dsTMatrix-method}
\alias{lu,dspMatrix-method}
\alias{lu,dsyMatrix-method}
\alias{lu,dtCMatrix-method}
\alias{lu,dtRMatrix-method}
\alias{lu,dtTMatrix-method}
\alias{lu,dtpMatrix-method}
\alias{lu,dtrMatrix-method}
\alias{lu,matrix-method}
\alias{lu,sparseMatrix-method}
%
\description{
  Computes the pivoted LU factorization of an \eqn{m \times n}{m-by-n}
  real matrix \eqn{A}, which has the general form
  \deqn{P_{1} A P_{2} = L U}{P1 * A * P2 = L * U}
  or (equivalently)
  \deqn{A = P_{1}' L U P_{2}'}{A = P1' * L * U * P2'}
  where
  \eqn{P_{1}}{P1} is an \eqn{m \times m}{m-by-m} permutation matrix, 
  \eqn{P_{2}}{P2} is an \eqn{n \times n}{n-by-n} permutation matrix, 
  \eqn{L} is an \eqn{m \times \text{min}(m,n)}{m-by-min(m,n)}
  unit lower trapezoidal matrix, and
  \eqn{U} is a  \eqn{\text{min}(m,n) \times n}{min(m,n)-by-n}
  upper trapezoidal matrix.

  Methods for \code{\linkS4class{denseMatrix}} are built on
  LAPACK routine \code{DGETRF}, which does not permute columns,
  so that \eqn{P_{2}}{P2} is an identity matrix.

  Methods for \code{\linkS4class{sparseMatrix}} are built on
  CSparse routine \code{cs_lu}, which requires \eqn{m = n},
  so that \eqn{L} and \eqn{U} are triangular matrices.
}
\usage{
lu(x, \dots)
\S4method{lu}{dgeMatrix}(x, warnSing = TRUE, \dots)
\S4method{lu}{dgCMatrix}(x, errSing = TRUE, order = NA_integer_, tol = 1,
   keep.dimnames = TRUE, \dots)
\S4method{lu}{dsyMatrix}(x, cache = TRUE, \dots)
\S4method{lu}{dsCMatrix}(x, cache = TRUE, \dots)
\S4method{lu}{matrix}(x, \dots)
}
\arguments{
  \item{x}{a \link[=is.finite]{finite} matrix or
    \code{\linkS4class{Matrix}} to be factorized,
    which must be square if sparse.}
  \item{warnSing}{a logical indicating if a \link{warning} should
    be signaled for singular \code{x}.  Used only by methods for
    dense matrices.}
  \item{errSing}{a logical indicating if an \link[=stop]{error} should
    be signaled for singular \code{x}.  Used only by methods for
    sparse matrices.}
  \item{order}{an integer in \code{0:3} passed to CSparse routine
    \code{cs_sqr}, indicating a strategy for choosing the column
    permutation \eqn{P_{2}}{P2}. 0 means no column permutation.
    1, 2, and 3 indicate a fill-reducing ordering of \eqn{A + A'},
    \eqn{\tilde{A}' \tilde{A}}{A~' * A~}, and \eqn{A' A}{A' * A},
    where \eqn{\tilde{A}}{A~} is \eqn{A} with \dQuote{dense} rows
    removed.
    \code{NA} (the default) is equivalent to 2 if \code{tol == 1}
    and 1 otherwise.
    Do not set to 0 unless you know that the column order of \eqn{A}
    is already sensible.}
  \item{tol}{a number.  The original pivot element is used
    if its absolute value exceeds \code{tol * a},
    where \code{a} is the maximum in absolute value of the
    other possible pivot elements.
    Set \code{tol < 1} only if you know what you are doing.}
  \item{keep.dimnames}{a logical indicating if \code{x@Dimnames}
    should be propagated to the result as \code{lu(x)@Dimnames}.}
  \item{cache}{a logical indicating if the result should be
    cached in \code{x@factors[["LU"]]}.  Note that
    caching is experimental and that only methods for classes
    extending \code{\linkS4class{compMatrix}} will have this
    argument.}
  \item{\dots}{further arguments passed to or from methods.}
}
\value{
  An object representing the factorization, inheriting from
  virtual class \code{\linkS4class{LU}}.  The specific class
  is \code{\linkS4class{denseLU}} unless \code{x} inherits
  from virtual class \code{\linkS4class{sparseMatrix}},
  in which case it is \code{\linkS4class{sparseLU}}.
}
\details{
  What happens when \code{x} is determined to be near-singular
  differs by method.  The method for class \code{\linkS4class{dgeMatrix}}
  completes the factorization, warning if \code{warnSing = TRUE}
  and in any case returning a valid \code{\linkS4class{denseLU}}
  object.  Users of this method can detect singular \code{x} with
  a suitable warning handler; see \code{\link{tryCatch}}.
  In contrast, the method for class \code{\linkS4class{dgCMatrix}}
  abandons further computation, throwing an error if \code{errSing = TRUE}
  and otherwise returning \code{NA}.  Users of this method can
  detect singular \code{x} with an error handler or by setting
  \code{errSing = FALSE} and testing for a formal result with
  \code{\link{isS4}}.

  Due to caching, it is currently unsafe to mix different values of
  \code{errSing}.  Notably, if \code{x} is a singular \code{dgCMatrix},
  then the second call below will return \code{NA} instead of
  throwing an error:
  \preformatted{    lu(x, errSing = FALSE); lu(x, errSing = TRUE)}
}
\seealso{
  Classes \code{\linkS4class{denseLU}} and
  \code{\linkS4class{sparseLU}} and their methods.

  Classes \code{\linkS4class{dgeMatrix}} and
  \code{\linkS4class{dgCMatrix}}.

  Generic function \code{\link{expand2}},
  for constructing matrix factors from the result.
  
  Generic functions \code{\link{Cholesky}}, \code{\link{BunchKaufman}},
  \code{\link{Schur}}, and \code{\link{qr}},
  for computing other factorizations.
}
\references{
  The LAPACK source code, including documentation; see
  \url{https://netlib.org/lapack/double/dgetrf.f}.
  
  Davis, T. A. (2006).
  \emph{Direct Methods for Sparse Linear Systems},
  Society for Industrial and Applied Mathematics.
}
\examples{
showMethods("lu", inherited = FALSE)

##--- Dense  -------------------------
x <- Matrix(rnorm(9), 3, 3)
lu(x)
dim(x2 <- round(10 * x[,-3]))# non-square
expand(lu2 <- lu(x2))

##--- Sparse (see more in ?"sparseLU-class")----- % ./sparseLU-class.Rd

pm <- as(readMM(system.file("external/pores_1.mtx",
                            package = "Matrix")),
         "CsparseMatrix")
str(pmLU <- lu(pm))		# p is a 0-based permutation of the rows
                                # q is a 0-based permutation of the columns
## permute rows and columns of original matrix
ppm <- pm[pmLU@p + 1L, pmLU@q + 1L]
pLU <- drop0(pmLU@L \%*\% pmLU@U) # L \%*\% U -- dropping extra zeros
## equal up to "rounding"
ppm[1:14, 1:5]
pLU[1:14, 1:5]
}
