\name{update-methods}
%%
%% \alias{update}
\alias{update-methods}
\alias{update,sparseCholesky-method}
%%
\title{Update a Sparse Cholesky Factorization}
\description{
  S4 methods defined in package \pkg{Matrix} for S4 generic function
  \code{update}, generalizing \code{stats::\link{update}}.

  The method for virtual class \code{\linkS4class{sparseCholesky}}
  adjusts the numeric entries of a sparse Cholesky factorization to a
  new sparse matrix.  It changes factorizations of \eqn{A}
  \deqn{P_1 A P_1' = L_1 D L_1' \overset{D_{jj} \ge 0}{=} L L'}{P1 * A * P1' = L1 * D * L1' [ = L * L' ]}
  to factorizations of \eqn{\tilde{A}}{A~}
  \deqn{P_1 \tilde{A} P_1' = \tilde{L}_1 \tilde{D} \tilde{L}_1' \overset{\tilde{D}_{jj} \ge 0}{=} \tilde{L} \tilde{L}'}{P1 * A~ * P1' = L~1 * D~ * L~1' [ = L~ * L~' ]}
  where the new matrix \eqn{\tilde{A}}{A~} is constrained to have
  nonzero pattern identical to that of the old matrix \eqn{A}.
}
\usage{
\S4method{update}{sparseCholesky}(object, parent, mult = 0, \dots)
}
\arguments{
  \item{object}{a \code{\linkS4class{sparseCholesky}} object.}
  \item{parent}{a \code{\linkS4class{sparseMatrix}} object, coerced to
    \code{\linkS4class{CsparseMatrix}} as necessary.}
  \item{mult}{a numeric or complex vector of positive length.}
  \item{\dots}{optional arguments not yet used.}
}
\value{
  A new \code{\linkS4class{sparseCholesky}} object representing the
  Cholesky factorization of \code{op(parent)} plus \code{mult[1]} times
  the identity matrix, where \code{op} is \code{\link{identity}} for
  formally Hermitian \code{parent} and \code{\link{tcrossprod}} for all
  other \code{parent}.
}
\seealso{
  Virtual class \code{\linkS4class{sparseCholesky}}.

  Generic functions \code{\linkS4methods{Cholesky}} and
  \code{\linkS4methods{updown}}.
}
\examples{
set.seed(44L)
X <- rsparsematrix(1000L, 200L, 0.01)
A <- crossprod(X)

u0 <- u0. <- Cholesky(A, LDL = FALSE)
u1 <- update(u0, parent = 4 * A)
u2 <- update(u0, parent =     A , mult = 2)
u3 <- update(u0, parent =   t(X), mult = 2)

u0.@x <- 2 * u0.@x

stopifnot(all.equal(u1, u0., tolerance = 1e-14),
          all.equal(u3, u2 , tolerance = 1e-14))
}
\docType{methods}
\keyword{algebra}
\keyword{array}
\keyword{methods}
