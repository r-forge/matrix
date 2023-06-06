\name{Schur-class}
\title{Schur Factorizations}
%
\docType{class}
\keyword{algebra}
\keyword{array}
\keyword{classes}
%
\alias{Schur-class}
%
\alias{determinant,Schur,logical-method}
%
\description{
  \code{Schur} is the class of Schur factorizations of
  \eqn{n \times n}{n-by-n} real matrices \eqn{A},
  having the general form
  \deqn{A = Q T Q'}{A = Q * T * Q'}
  where
  \eqn{Q} is an orthogonal matrix and
  \eqn{T} is a block upper triangular matrix with
  \eqn{1 \times 1}{1-by-1} or \eqn{2 \times 2}{2-by-2} diagonal blocks
  specifying the real and complex conjugate eigenvalues of \eqn{A}.
  The column vectors of \eqn{Q} are the Schur vectors of \eqn{A},
  and \eqn{T} is the Schur form of \eqn{A}.

  The Schur factorization generalizes the spectral decomposition
  of normal matrices \eqn{A}, whose Schur form is block diagonal,
  to arbitrary square matrices.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}}{inherited from virtual class
      \code{\linkS4class{MatrixFactorization}}.}
    \item{\code{Q}}{an orthogonal matrix,
      inheriting from virtual class \code{\linkS4class{Matrix}}.}
    \item{\code{T}}{a block upper triangular matrix,
      inheriting from virtual class \code{\linkS4class{Matrix}}.
      The diagonal blocks have dimensions 1-by-1 or 2-by-2.}
    \item{\code{EValues}}{a numeric or complex vector containing
      the eigenvalues of the diagonal blocks of \code{T}, which are
      the eigenvalues of \code{T} and consequently of the factorized
      matrix.}
  }
}
\section{Extends}{
  Class \code{\linkS4class{SchurFactorization}}, directly.
  Class \code{\linkS4class{MatrixFactorization}}, by class
  \code{\linkS4class{SchurFactorization}}, distance 2.
}
\section{Instantiation}{
  Objects can be generated directly by calls of the form
  \code{new("Schur", ...)}, but they are more typically obtained
  as the value of \code{\link{Schur}(x)} for \code{x} inheriting from
  \code{\linkS4class{Matrix}} (often \code{\linkS4class{dgeMatrix}}).
}
\section{Methods}{
  \describe{
    \item{\code{determinant}}{\code{signature(from = "Schur", logarithm = "logical")}:
      computes the determinant of the factorized matrix \eqn{A}
      or its logarithm.}
    \item{expand2}{\code{signature(x = "Schur")}:
      see \code{\link{expand2-methods}}.}
    \item{\code{solve}}{\code{signature(a = "Schur", b = .)}:
      see \code{\link{solve-methods}}.}
  }
}
\details{
  The matrix \eqn{A} and its Schur form \eqn{T} are \emph{similar}
  and thus have the same spectrum.  The eigenvalues are computed
  trivially as the eigenvalues of the diagonal blocks of \eqn{T}.
}
\seealso{
  Class \code{\linkS4class{dgeMatrix}}.
  
  Generic functions \code{\link{expand2}} and \code{\link{Schur}}.
}
\references{
  The LAPACK source code, including documentation; see
  \url{https://netlib.org/lapack/double/dgees.f}.
}
\examples{
showClass("Schur")
set.seed(0)

(A <- Matrix(rnorm(16L), 4L, 4L))
(Sch.A <- Schur(A))
(e.Sch.A <- expand2(Sch.A))

## A ~ Q T Q' in floating point
stopifnot(all.equal(A, with(e.Sch.A, Q \%*\% T \%*\% Q.)))

## Factorization handled as factorized matrix
b <- rnorm(4L)
stopifnot(all.equal(det(A), det(Sch.A)),
          all.equal(solve(A, b), solve(Sch.A, b)))

## One of the non-general cases:
Schur(Diagonal(5L))
}
