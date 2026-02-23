\name{updown-methods}
%%
\alias{updown}
\alias{updown-methods}
\alias{updown,character,ANY,ANY-method}
\alias{updown,logical,Matrix,sparseCholesky-method}
\alias{updown,logical,matrix,sparseCholesky-method}
%%
\title{Rank-k Updates and Downdates of Cholesky Factorizations}
\description{
  Computes a rank-\eqn{k} update or downdate of a Cholesky factorization
  \deqn{P_{1} A P_{1}' = L_{1} D L_{1}' = L L'}{P1 * A * P1' = L1 * D * L1' = L * L'}
  which for some \eqn{k}-column matrix \eqn{C} is the factorization
  \deqn{P_{1} (A + s C C') P_{1}' = \tilde{L}_{1} \tilde{D} \tilde{L}_{1}' = \tilde{L} \tilde{L}'}{P1 * (A + s * C * C') * P1' = L~1 * D~ * L~1' = L~ * L~'}
  Here, \eqn{s = 1} for an update and \eqn{s = -1} for a downdate.
}
\usage{
updown(update, C, L)

\S4method{updown}{character,ANY,ANY}(update, C, L)
\S4method{updown}{logical,Matrix,sparseCholesky}(update, C, L)
\S4method{updown}{logical,matrix,sparseCholesky}(update, C, L)
}
\arguments{
  \item{update}{a logical (\code{TRUE} or \code{FALSE}) or character
    (\code{"+"} or \code{"-"}) indicating if \code{L} should be updated
    (or otherwise downdated).}
  \item{C}{an \R{} object representing the \eqn{k}-column matrix
    \eqn{C}, such that \code{\link{tcrossprod}(C)} has the dimensions of
    \code{L}.  Package \pkg{Matrix} provides methods for
    \code{\linkS4class{Matrix}} objects and traditional matrices;
    these coerce to \code{\linkS4class{CsparseMatrix}} as necessary.}
  \item{L}{an \R{} object representing the original Cholesky
    factorization.  Package \pkg{Matrix} provides methods for virtual
    class \code{\linkS4class{sparseCholesky}}.}
}
\value{
  An \R{} object representing the updated or downdated Cholesky
  factorization.

  Methods provided by package \pkg{Matrix} return a new
  \code{\linkS4class{sparseCholesky}} object.  The new factorization is
  square root-free, even if the old one is not, hence it is required
  that \eqn{A + s C C'}{A + s * C * C'} has nonzero leading principal
  minors, not that it is positive definite.
}
\seealso{
  Virtual class \code{\linkS4class{sparseCholesky}}.

  Generic function \code{\linkS4methods{Cholesky}}.

  The method for generic function \code{\linkS4methods{update}} with
  signature \code{c(object="sparseCholesky")}, which is \emph{not}
  equivalent to \code{updown(update = TRUE, C = ., L = object)}.
}
\author{
  Nicholas Nagle (University of Tennessee, Knoxville) provided the
  original implementation.
}
\references{
  Davis, T. A., Hager, W. W. (2001).
  Multiple-rank modifications of a sparse Cholesky factorization.
  \emph{SIAM Journal on Matrix Analysis and Applications},
  \emph{22}(4), 997-1013.
  \doi{10.1137/S0895479899357346}
}
\examples{
X  <- new("dgCMatrix", Dim = c(3L, 5L), p = c(0L, 1L, 2L, 4L, 4L, 6L),
          i = c(2L, 0L, 1L, 2L, 0L, 1L), x = as.double(1:6))
I5 <- new("ddiMatrix", Dim = c(5L, 5L), diag = "U")
A <- crossprod(X) + I5
u0 <- Cholesky(A)
u1 <- updown("+", I5, u0)
u2 <- updown("-", I5, u1)
stopifnot(all.equal(u0, u2))
}
\docType{methods}
\keyword{algebra}
\keyword{array}
\keyword{methods}
