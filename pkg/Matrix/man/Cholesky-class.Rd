\name{Cholesky-class}
\title{Dense Cholesky Factorizations}
%
\docType{class}
\keyword{algebra}
\keyword{array}
\keyword{classes}
%
\alias{Cholesky-class}
\alias{pCholesky-class}
%
\alias{coerce,Cholesky,dtrMatrix-method}
\alias{coerce,pCholesky,dtpMatrix-method}
\alias{determinant,Cholesky,logical-method}
\alias{determinant,pCholesky,logical-method}
\alias{diag,Cholesky-method}
\alias{diag,pCholesky-method}
%
\description{
  Classes \code{Cholesky} and \code{pCholesky} represent
  dense Cholesky factorizations of \eqn{n \times n}{n-by-n}
  real, symmetric, positive semidefinite matrices \eqn{A},
  having the general form
  \deqn{P_{1} A P_{1}' = L_{1} D L_{1}' = L L'}{P1 * A * P1' = L1 * D * L1' = L * L'}
  or (equivalently)
  \deqn{A = P_{1}' L_{1} D L_{1}' P_{1} = P_{1}' L L' P_{1}}{A = P1' * L1 * D * L1' * P1 = P1' * L * L' * P1}
  where
  \eqn{P_{1}}{P1} is a permutation matrix,
  \eqn{L_{1}}{L1} is a unit lower triangular matrix,
  \eqn{D} is a non-negative diagonal matrix, and
  \eqn{L = L_{1} \sqrt{D}}{L = L1 * sqrt(D)}.

  These classes store the entries of the Cholesky factor
  \eqn{L} or its transpose \eqn{L'} in a dense format as
  a vector of length \eqn{nn}{n*n} (\code{Cholesky}) or
  \eqn{n(n+1)/2}{n*(n+1)/2} (\code{pCholesky}), the latter
  giving the \dQuote{packed} representation.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}}{inherited from virtual class
      \code{\linkS4class{MatrixFactorization}}.}
    \item{\code{uplo}}{a string, either \code{"U"} or \code{"L"},
      indicating which triangle (upper or lower) of the factorized
      symmetric matrix was used to compute the factorization and
      in turn whether \code{x} stores \eqn{L'} or \eqn{L}.}
    \item{\code{x}}{a numeric vector of length \code{n*n}
      (\code{Cholesky}) or \code{n*(n+1)/2} (\code{pCholesky}),
      where \code{n=Dim[1]}, listing the entries of the Cholesky
      factor \eqn{L} or its transpose \eqn{L'} in column-major
      order.}
    \item{\code{perm}}{a 1-based integer vector of length \code{Dim[1]}
      specifying the permutation applied to the rows and columns
      of the factorized matrix.  \code{perm} of length 0 is valid and
      equivalent to the identity permutation, implying no pivoting.}
  }
}
\section{Extends}{
  Class \code{\linkS4class{CholeskyFactorization}}, directly.
  Class \code{\linkS4class{MatrixFactorization}}, by class
  \code{\linkS4class{CholeskyFactorization}}, distance 2.
}
\section{Instantiation}{
  Objects can be generated directly by calls of the form
  \code{new("Cholesky", ...)} or \code{new("pCholesky", ...)},
  but they are more typically obtained as the value of
  \code{\link{Cholesky}(x)} for \code{x} inheriting from
  \code{\linkS4class{dsyMatrix}} or \code{\linkS4class{dspMatrix}}
  (often the subclasses of those reserved for positive
  semidefinite matrices, namely \code{\linkS4class{dpoMatrix}}
  and \code{\linkS4class{dppMatrix}}).
}
\section{Methods}{
  \describe{
    \item{\code{coerce}}{\code{signature(from = "Cholesky", to = "dtrMatrix")}:
      returns a \code{\linkS4class{dtrMatrix}} representing
      the Cholesky factor \eqn{L} or its transpose \eqn{L'};
      see \sQuote{Note}.}
    \item{\code{coerce}}{\code{signature(from = "pCholesky", to = "dtpMatrix")}:
      returns a \code{\linkS4class{dtpMatrix}} representing
      the Cholesky factor \eqn{L} or its transpose \eqn{L'};
      see \sQuote{Note}.}
    \item{\code{determinant}}{\code{signature(from = "p?Cholesky", logarithm = "logical")}:
      computes the determinant of the factorized matrix \eqn{A}
      or its logarithm.}
    \item{\code{diag}}{\code{signature(x = "p?Cholesky")}:
      returns a numeric vector of length \eqn{n} containing the diagonal
      elements of \eqn{D}, which are the squared diagonal elements of
      \eqn{L}.}
    \item{\code{expand2}}{\code{signature(x = "p?Cholesky")}:
      see \code{\link{expand2-methods}}.}
    \item{\code{solve}}{\code{signature(a = "p?Cholesky", b = .)}:
      see \code{\link{solve-methods}}.}
  }
}
\note{
  In \pkg{Matrix} \code{< 1.6-0}, class \code{Cholesky} extended
  \code{\linkS4class{dtrMatrix}} and class \code{pCholesky} extended
  \code{\linkS4class{dtpMatrix}}, reflecting the fact that the factor
  \eqn{L} is indeed a triangular matrix.
  \pkg{Matrix} \code{1.6-0} removed these extensions so that methods
  would no longer be inherited from \code{dtrMatrix} and \code{dtpMatrix}.  
  The availability of such methods gave the wrong impression that
  \code{Cholesky} and \code{pCholesky} represent proper matrices,
  when in fact they represent an ordered set of matrix factors.
  
  The coercions \code{as(., "dtrMatrix")} and \code{as(., "dtpMatrix")}
  are provided for users who understand the caveats.
}
\seealso{
  Class \code{\linkS4class{CHMfactor}} for sparse Cholesky factorizations.

  Classes \code{\linkS4class{dpoMatrix}} and \code{\linkS4class{dppMatrix}}.
  
  Generic functions \code{\link{expand2}} and \code{\link{Cholesky}}.
}
\examples{
(sp <- pack(Matrix(diag(5L) + 1)))
signif(chsp <- chol(sp), 4L)

(sy <- crossprod(Matrix(rnorm(18L), 6L, 3L)))
signif(chsy <- chol(sy), 4L)

(sy. <- if (chsy@uplo == "U") crossprod(chsy) else tcrossprod(chsy))
stopifnot(all.equal(as(sy., "matrix"), as(sy, "matrix"),
                    tolerance = 1e-14))
}
