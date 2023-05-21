\name{qr-methods}
\title{Methods for QR Factorization}
%
\docType{methods}
\keyword{algebra}
\keyword{array}
\keyword{methods}
%
\alias{qr}
\alias{qr-methods}
%
\alias{qr,dgCMatrix-method}
\alias{qr,sparseMatrix-method}
%
\description{
  Computes the pivoted QR factorization of an \eqn{m \times n}{m-by-n}
  real matrix \eqn{A}, which has the general form
  \deqn{P_{1} A P_{2} = Q R}{P1 * A * P2 = Q * R}
  or (equivalently)
  \deqn{A = P_{1}' Q R P_{2}'}{A = P1' * Q * R * P2'}
  where
  \eqn{P_{1}}{P1} and \eqn{P_{2}}{P2} are permutation matrices,
  \eqn{Q} is an \eqn{m \times m}{m-by-m} orthogonal matrix, and
  \eqn{R} is an \eqn{m \times n}{m-by-n} upper trapezoidal matrix.

  \code{\linkS4class{denseMatrix}} use the default method implemented
  in \pkg{base}, namely \code{\link{qr.default}}.  It is built on
  LINPACK routine \code{dqrdc} and LAPACK routine \code{dgeqp3}, which
  do not pivot rows, so that \eqn{P_{1}}{P1} is an identity matrix.
  
  Methods for \code{\linkS4class{sparseMatrix}} are built on
  CSparse routines \code{cs_sqr} and \code{cs_qr}, which require
  \eqn{m \ge n}{m >= n}.
}
\usage{
qr(x, \dots)
\S4method{qr}{dgCMatrix}(x, order = 3L, \dots)
}
\arguments{
  \item{x}{a \link[=is.finite]{finite} matrix or
    \code{\linkS4class{Matrix}} to be factorized,
    satisfying \code{nrow(x) >= ncol(x)} if sparse.}
  \item{order}{an integer in \code{0:3} passed to CSparse routine
    \code{cs_sqr}, indicating a strategy for choosing the column
    permutation \eqn{P_{2}}{P2}.  0 means no column permutation.
    1, 2, and 3 indicate a fill-reducing ordering of \eqn{A + A'},
    \eqn{\tilde{A}' \tilde{A}}{A~' * A~}, and \eqn{A' A}{A' * A},
    where \eqn{\tilde{A}}{A~} is \eqn{A} with \dQuote{dense} rows
    removed.
    Do not set to 0 unless you know that the column order of \eqn{A}
    is already sensible.}
  \item{\dots}{further arguments passed to or from methods.}
}
\value{
  An object representing the factorization, of S4 class
  \code{\linkS4class{sparseQR}} or S3 class \code{\link[base]{qr}}.
  \code{sparseQR} is given for \code{x} inheriting from virtual class
  \code{\linkS4class{sparseMatrix}}.  \code{qr} is given for \code{x}
  inheriting from virtual class \code{\linkS4class{denseMatrix}},
  as well as for traditional matrices.
}
\details{
  If \code{x} is sparse and structurally rank deficient, having
  structural rank \eqn{r < n}, then \code{x} is augmented with
  \eqn{(n-r)} rows of (partly non-structural) zeros, such that
  the augmented matrix has structural rank \eqn{n}.
  This augmented matrix is factorized as described above:
  \deqn{P_{1} A P_{2} = P_{1} \begin{bmatrix} A_{0} \\ 0 \end{bmatrix} P_{2} = Q R}{P1 * A * P2 = P1 * [A0; 0] * P2 = Q * R}
  where \eqn{A_{0}} denotes the original, user-supplied
  \eqn{(m-(n-r)) \times n}{(m-(n-r))-by-n} matrix.
}
\seealso{
  Class \code{\linkS4class{sparseQR}} and its methods.

  Class \code{\linkS4class{dgCMatrix}}.
  
  Generic function \code{\link[base]{qr}} from \pkg{base},
  whose default method \code{qr.default} \dQuote{defines}
  the S3 class \code{qr} of dense QR factorizations.

  Generic function \code{\link{expand2}},
  for constructing matrix factors from the result.
  
  Generic functions \code{\link{Cholesky}}, \code{\link{BunchKaufman}},
  \code{\link{Schur}}, and \code{\link{lu}},
  for computing other factorizations.
}
\references{
  Davis, T. A. (2006).
  \emph{Direct Methods for Sparse Linear Systems},
  Society for Industrial and Applied Mathematics.
}
\examples{
showMethods("qr", inherited = FALSE)

%% FIXME: Currently mixing example + regression tests <--> ../tests/factorizing.R
##------------- example of pivoting -- from base'  qraux.Rd -------------
X <- cbind(int = 1,
           b1=rep(1:0, each=3), b2=rep(0:1, each=3),
           c1=rep(c(1,0,0), 2), c2=rep(c(0,1,0), 2), c3=rep(c(0,0,1),2))
rownames(X) <- paste0("r", seq_len(nrow(X)))
dnX <- dimnames(X)
bX <- X # [b]ase version of X
X <- as(bX, "sparseMatrix")
X # is singular, columns "b2" and "c3" are "extra"
stopifnot(identical(dimnames(X), dnX))# some versions changed X's dimnames!
c(rankMatrix(X)) # = 4 (not 6)
m <- function(.) as(., "matrix")

##----- regular case ------------------------------------------
Xr <- X[ , -c(3,6)] # the "regular" (non-singular) version of X
stopifnot(rankMatrix(Xr) == ncol(Xr))
Y <- cbind(y <- setNames(1:6, paste0("y", 1:6)))
## regular case:
qXr   <- qr(  Xr)
qxr   <- qr(m(Xr))
qxrLA <- qr(m(Xr), LAPACK=TRUE) # => qr.fitted(), qr.resid() not supported
qcfXy <- qr.coef (qXr, y) # vector
qcfXY <- qr.coef (qXr, Y) # 4x1 dgeMatrix
cf <- c(int=6, b1=-3, c1=-2, c2=-1)
doExtras <- interactive() || nzchar(Sys.getenv("R_MATRIX_CHECK_EXTRA"))
tolE <- if(doExtras) 2e-15 else 1e-13
stopifnot(exprs = {
  all.equal(qr.coef(qxr,  y),   cf,  tol=tolE)
  all.equal(qr.coef(qxrLA,y),   cf,  tol=tolE)
  all.equal(qr.coef(qxr,  Y), m(cf), tol=tolE)
  all.equal(  qcfXy,    cf, tol=tolE)
  all.equal(m(qcfXY), m(cf), tol=tolE)
  all.equal(y, qr.fitted(qxr, y), tol=2*tolE)
  all.equal(y, qr.fitted(qXr, y), tol=2*tolE)
  all.equal(m(qr.fitted(qXr, Y)), qr.fitted(qxr, Y), tol=tolE)
  all.equal(  qr.resid (qXr, y),  qr.resid (qxr, y), tol=tolE)
  all.equal(m(qr.resid (qXr, Y)), qr.resid (qxr, Y), tol=tolE)
})

##----- rank-deficient ("singular") case ------------------------------------

(qX <- qr(  X))           # both @p and @q are non-trivial permutations
 qx <- qr(m(X)) ; str(qx) # $pivot is non-trivial, too

drop0(R. <- qr.R(qX), tol=tolE) # columns *permuted*: c3 b1 ..
Q. <- qr.Q(qX)
qI <- sort.list(qX@q) # the inverse 'q' permutation
(X. <- drop0(Q. \%*\% R.[, qI], tol=tolE))## just = X, incl. correct colnames
stopifnot(all(X - X. < 8*.Machine$double.eps),
          ## qrR(.) returns R already "back permuted" (as with qI):
          identical(R.[, qI], qrR(qX, row.names = FALSE)) )
##
## In this sense, classical qr.coef() is fine:
cfqx <- qr.coef(qx, y) # quite different from
nna <- !is.na(cfqx)
stopifnot(all.equal(unname(qr.fitted(qx,y)),
                    as.numeric(X[,nna] \%*\% cfqx[nna])))
## FIXME: do these make *any* sense? --- should give warnings !
qr.coef(qX, y)
qr.coef(qX, Y)
rm(m)
}
