\name{rcond-methods}
\title{Estimate the (Reciprocal) Condition Number}
%
\docType{methods}
\keyword{algebra}
\keyword{math}
\keyword{methods}
%
\alias{kappa}
\alias{kappa-methods}
\alias{rcond}
\alias{rcond-methods}
%
\alias{kappa,Matrix-method}
%
\alias{rcond,ANY,missing-method}
\alias{rcond,denseMatrix,character-method}
\alias{rcond,diagonalMatrix,character-method}
\alias{rcond,indMatrix,character-method}
\alias{rcond,pMatrix,character-method}
\alias{rcond,sparseMatrix,character-method}
%
\usage{
kappa(z, \dots)
rcond(x, norm, \dots)

\S4method{rcond}{sparseMatrix,character}(x, norm, exact = FALSE, inverse = solve(x),
      warn = TRUE, \dots)
}
\description{
  Estimate the reciprocal of the condition number of a matrix.

  This is a generic function with several methods, as seen by
  \code{\link{showMethods}(rcond)}.
}
\arguments{
  \item{x, z}{an \R object that inherits from the \code{Matrix} class.}
  \item{norm}{character string indicating the type of norm to be used in
    the estimate. The default is \code{"O"} for the 1-norm (\code{"O"} is
    equivalent to \code{"1"}).  For sparse matrices, when \code{useInv=TRUE},
    \code{norm} can be any of the \code{kind}s allowed for \code{\link{norm}};
    otherwise, the other possible value is \code{"I"} for the infinity
    norm, see also \code{\link{norm}}.
  }
  \item{exact}{.}
  \item{inverse}{.}
  \item{warn}{.}
  %% \item{useInv}{logical (or \code{"Matrix"} containing
  %%   \code{\link{solve}(x)}).  If not false, compute the reciprocal
  %%   condition number as \eqn{1/(\|x\| \cdot \|x^{-1}\|)}{1/(||x|| * ||x^(-1)||)},
  %%   where \eqn{x^{-1}}{x^(-1)} is the inverse of \eqn{x}, \code{solve(x)}.

  %%   This may be an efficient alternative (only) in situations where
  %%   \code{solve(x)} is fast (or known), e.g., for (very) sparse or
  %%   triangular matrices.

  %%   Note that the \emph{result} may differ depending on \code{useInv},
  %%   as per default, when it is false, an \emph{approximation} is
  %%   computed.
  %% }
  \item{\dots}{further arguments passed to or from other methods.}
}
\value{
  An estimate of the reciprocal condition number of \code{x}.
}
\section{BACKGROUND}{
  The condition number of a regular (square) matrix is the product of
  the \code{\link{norm}} of the matrix and the norm of its inverse (or
  pseudo-inverse).

  More generally, the condition number is defined (also for
  non-square matrices \eqn{A}) as
  \deqn{\kappa(A) = \frac{\max_{\|v\| = 1} \|A v\|}{\min_{\|v\| = 1} \|A v\|}.}{%
    kappa(A) = (max_(||v|| = 1; || Av ||)) /(min_(||v|| = 1; || Av ||)).}
  Whenever \code{x} is \emph{not} a square matrix, in our method
  definitions, this is typically computed via \code{rcond(qr.R(qr(X)), ...)}
  where \code{X} is \code{x} or \code{t(x)}.

  The condition number takes on values between 1 and infinity,
  inclusive, and can be viewed as a factor by which errors in solving
  linear systems with this matrix as coefficient matrix could be
  magnified.

  \code{rcond()} computes the \emph{reciprocal} condition number
  \eqn{1/\kappa} with values in \eqn{[0,1]} and can be viewed as a
  scaled measure of how close a matrix is to being rank deficient (aka
  \dQuote{singular}).

  Condition numbers are usually estimated, since exact computation is
  costly in terms of floating-point operations.  An (over) estimate of
  reciprocal condition number is given, since by doing so overflow is
  avoided.  Matrices are well-conditioned if the reciprocal condition
  number is near 1 and ill-conditioned if it is near zero.
}
\seealso{
  \code{\link{norm}}, \code{\link[base]{kappa}()} from package
  \pkg{base} computes an \emph{approximate} condition number of a
  \dQuote{traditional} matrix, even non-square ones, with respect to the
  \eqn{p=2} (Euclidean) \code{\link{norm}}.
  \code{\link[base]{solve}}.

  \code{\link{condest}}, a newer \emph{approximate} estimate of
  the (1-norm) condition number, particularly efficient for large sparse
  matrices.
}
\references{
  Golub, G., and Van Loan, C. F. (1989).
  \emph{Matrix Computations,}
  2nd edition, Johns Hopkins, Baltimore.
}
\examples{
\dontshow{ % for R_DEFAULT_PACKAGES=NULL
library(stats, pos = "package:base", verbose = FALSE)
}
x <- Matrix(rnorm(9), 3, 3)
rcond(x)
## typically "the same" (with more computational effort):
1 / (norm(x) * norm(solve(x)))
rcond(Hilbert(9))  # should be about 9.1e-13

## For non-square matrices:
rcond(x1 <- cbind(1,1:10))# 0.05278
rcond(x2 <- cbind(x1, 2:11))# practically 0, since x2 does not have full rank

## sparse
(S1 <- Matrix(rbind(0:1,0, diag(3:-2))))
rcond(S1)
m1 <- as(S1, "denseMatrix")
all.equal(rcond(S1), rcond(m1))

## wide and sparse
rcond(Matrix(cbind(0, diag(2:-1))))

## Large sparse example ----------
m <- Matrix(c(3,0:2), 2,2)
M <- bdiag(kronecker(Diagonal(2), m), kronecker(m,m))
36*(iM <- solve(M)) # still sparse
MM <- kronecker(Diagonal(10), kronecker(Diagonal(5),kronecker(m,M)))
dim(M3 <- kronecker(bdiag(M,M),MM)) # 12'800 ^ 2
if(interactive()) ## takes about 2 seconds if you have >= 8 GB RAM
  system.time(r <- rcond(M3))
## whereas this is *fast* even though it computes  solve(M3)
system.time(r. <- rcond(M3, exact=TRUE))
if(interactive()) ## the values are not the same
  c(r, r.)  # 0.05555 0.013888
## for all norms available for sparseMatrix, except "2" :
cbind(rr <- sapply(c("O","I","F","M"),
             function(N) rcond(M3, norm=N, exact=TRUE)))
\dontshow{stopifnot(all.equal(r., 1/72, tolerance=1e-12))}
}

%% \item{rcond}{\code{signature(x = "dpoMatrix", norm = "character")}:
%%   Returns (and stores) the reciprocal of the condition number of
%%   \code{x}.  The \code{norm} can be \code{"O"} for the
%%   one-norm (the default) or \code{"I"} for the infinity-norm.  For
%%   symmetric matrices the result does not depend on the norm.}
