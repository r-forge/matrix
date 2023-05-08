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
\alias{determinant,Cholesky,logical-method}
\alias{determinant,pCholesky,logical-method}
\alias{show,Cholesky-method}
\alias{show,pCholesky-method}
%
\description{
  Classes \code{Cholesky} and \code{pCholesky} represent
  Cholesky factorizations of \eqn{n \times n}{n-by-n} real,
  symmetric, positive semidefinite matrices \eqn{A},
  having the general form
  \deqn{A = L_{1} D L_{1}' = L L'}{A = L1 * D * L1' = L * L'}
  where
  \eqn{L_{1}}{L1} is a unit lower triangular matrix,
  \eqn{D} is a non-negative diagonal matrix, and
  \eqn{L = L_{1} \sqrt{D}}{L = L1 * sqrt(D)}.

  These classes store the entries of the Cholesky factor \eqn{L}
  in a dense format as a vector of length
  \eqn{nn}{n*n} (\code{Cholesky}) or
  \eqn{n(n+1)/2}{n*(n+1)/2} (\code{pCholesky}),
  the latter giving the \dQuote{packed} representation.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}}{inherited from virtual class
      \code{\linkS4class{MatrixFactorization}}.}
    \item{\code{uplo}, \code{diag}}{inherited from virtual class
      \code{\linkS4class{triangularMatrix}}.}
    \item{\code{x}}{inherited from virtual class
      \code{\linkS4class{dMatrix}}.  It contains the entries of the
      upper or lower triangular Cholesky factor in column-major order.}
  }
}
\section{Extends}{
  Class \code{\linkS4class{CholeskyFactorization}}, directly.
  Class \code{\linkS4class{MatrixFactorization}}, by class
  \code{\linkS4class{CholeskyFactorization}}, distance 2.

  Class \code{\linkS4class{dtrMatrix}} or
  \code{\linkS4class{dtpMatrix}}, directly.
  Class \code{\linkS4class{Matrix}}, by class
  \code{\linkS4class{dtrMatrix}} or class
  \code{\linkS4class{dtpMatrix}}, distance 3.
}
\section{Instantiation}{
  Objects can be generated directly by calls of the form
  \code{new("Cholesky", ...)} or \code{new("pCholesky", ...)},
  but they are more typically obtained as the value of
  \code{\link{chol}(x)} for \code{x} inheriting from
  \code{\linkS4class{dsyMatrix}} or \code{\linkS4class{dspMatrix}}
  (often the subclasses of those reserved for positive
  semidefinite matrices, namely \code{\linkS4class{dpoMatrix}}
  and \code{\linkS4class{dppMatrix}}).
}
\section{Methods}{
  \describe{
    \item{\code{determinant}}{\code{signature(from = "p?Cholesky", logarithm = "logical")}:
      computes the determinant of the factorized matrix \eqn{A}
      or its logarithm.}
    \item{\code{expand2}}{\code{signature(x = "p?Cholesky")}:
      see \code{\link{expand2-methods}}.}
    \item{\code{show}}{\code{signature(object = "p?Cholesky")}:
      prints the Cholesky factor as a typical
      \code{\linkS4class{Matrix}} rather than as a typical
      \code{\linkS4class{MatrixFactorization}},
      i.e., as a rectangular array with \code{Dimnames} along the margin,
      rather than as a list of slot values.}
  }
}
\note{
  For traditional numeric matrices \code{x}, \code{chol(x)} is a
  traditional numeric matrix and also \emph{triangular}.  Accordingly,
  classes \code{Cholesky} and \code{pCholesky} extend the subclasses
  \code{\linkS4class{dtrMatrix}} and \code{\linkS4class{dtpMatrix}}
  of virtual class \code{\linkS4class{denseMatrix}} and inherit almost
  all of the corresponding methods (for arithmetic, comparison, etc.).
  A corollary is that \code{determinant(R <- chol(A))} returns the
  determinant of \code{R}, not of \code{A}.
  
  Class \code{Cholesky} should not be confused with the
  \link[=Cholesky]{function of the same name}.  The latter returns
  objects inheriting from virtual class \code{\linkS4class{CHMfactor}}
  (rather than from class \code{Cholesky}) that specify the Cholesky
  factor in \emph{sparse} (rather than dense) format.
}
\seealso{
  Class \code{\linkS4class{CHMfactor}} for sparse Cholesky
  factorizations.

  Classes \code{\linkS4class{dpoMatrix}} and
  \code{\linkS4class{dtrMatrix}} and their packed counterparts.

  Generic functions \code{\link{expand2}}, \code{\link{chol}},
  and \code{\link{Cholesky}}.
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
