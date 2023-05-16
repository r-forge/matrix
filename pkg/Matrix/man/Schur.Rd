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
}
\examples{
Schur(Hilbert(9))              # Schur factorization (real eigenvalues)

(A <- Matrix(round(rnorm(5*5, sd = 100)), nrow = 5))
(Sch.A <- Schur(A))

eTA <- eigen(Sch.A@T)
str(SchA <- Schur(A, vectors=FALSE))# no 'T' ==> simple list
stopifnot(all.equal(eTA$values, eigen(A)$values, tolerance = 1e-13),
          all.equal(eTA$values,
                    local({z <- Sch.A@EValues
                           z[order(Mod(z), decreasing=TRUE)]}), tolerance = 1e-13),
          identical(SchA$T, Sch.A@T),
          identical(SchA$EValues, Sch.A@EValues))

## For the faint of heart, we provide Schur() also for traditional matrices:

a.m <- function(M) unname(as(M, "matrix"))
a <- a.m(A)
Sch.a <- Schur(a)
stopifnot(identical(Sch.a, list(Q = a.m(Sch.A @ Q),
				T = a.m(Sch.A @ T),
				EValues = Sch.A@EValues)),
	  all.equal(a, with(Sch.a, Q \%*\% T \%*\% t(Q)))
)
}
