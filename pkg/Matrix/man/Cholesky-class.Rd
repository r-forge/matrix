\name{Cholesky-class}
\title{Cholesky Factorizations}
\docType{class}
%
\alias{Cholesky-class}
\alias{pCholesky-class}
%
\alias{show,Cholesky-method}
\alias{show,pCholesky-method}
%
\description{
  Classes \code{Cholesky} and \code{pCholesky} represent Cholesky
  factorizations of real, symmetric, positive semidefinite matrices.
  They store the entries of the \code{n}-by-\code{n} Cholesky factor
  in dense format as a vector of length \code{n*n} (\code{Cholesky})
  or \code{n*(n+1)/2} (\code{pCholesky}), the latter giving the
  \dQuote{\bold{p}acked} representation.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}}{inherited from virtual class
      \code{\linkS4class{Matrix}}.}
    \item{\code{uplo}, \code{diag}}{inherited from virtual class
      \code{\linkS4class{triangularMatrix}}.}
    \item{\code{x}}{inherited from virtual class
      \code{\linkS4class{dMatrix}} and contains the entries of the
      upper or lower triangular Cholesky factor in column-major order.}
  }
}
\section{Extends}{
  Class \code{\linkS4class{CholeskyFactorization}}, directly.
  Class \code{\linkS4class{MatrixFactorization}}, by class
  \code{\linkS4class{CholeskyFactorization}}, distance 2.

  Class \code{\linkS4class{dtrMatrix}} or
  \code{\linkS4class{dtpMatrix}}, directly.
  Class \code{\linkS4class{denseMatrix}}, by class
  \code{\linkS4class{dtrMatrix}} or
  \code{\linkS4class{dtpMatrix}}, distance 3.
}
\section{Methods}{
  For traditional numeric matrices \code{x}, \code{chol(x)} is a
  traditional numeric matrix and also \emph{triangular}.  Accordingly,
  classes \code{Cholesky} and \code{pCholesky} extend the subclasses
  \code{\linkS4class{dtrMatrix}} and \code{\linkS4class{dtpMatrix}}
  of virtual class \code{\linkS4class{denseMatrix}} and inherit almost
  all of the corresponding methods (for arithmetic, comparison, etc.).

  A corollary is that \code{\link{determinant}(R <- chol(A))} returns
  the determinant of \code{R}, not of \code{A}.  This is currently the
  only exception to the rule that \code{\link{determinant}(x)} should
  return the determinant of the original, factorized matrix if \code{x}
  inherits from \code{\linkS4class{MatrixFactorization}}.
  
  Classes \code{Cholesky} and \code{pCholesky} are also unique within
  the virtual class \code{\linkS4class{MatrixFactorization}} in that
  they are \link{show}n as a typical \code{\linkS4class{Matrix}} would
  be, i.e., as a rectangular array rather than as a list of slot values.
}
\section{Objects from the class}{
  Objects can be created directly by calls of the form
  \code{new("Cholesky", ...)} or \code{new("pCholesky", ...)}.
  More typically, they are obtained as the value of
  \code{\link{chol}(x)} for \code{x} inheriting from
  \code{\linkS4class{dsyMatrix}} or \code{\linkS4class{dspMatrix}}
  (often the subclasses of those reserved for positive semidefinite
  matrices, namely \code{\linkS4class{dpoMatrix}} and
  \code{\linkS4class{dppMatrix}}).
}
\note{
  Class \code{Cholesky} should not be confused with the
  \link[=Cholesky]{function of the same name}.  The latter returns
  objects inheriting from virtual class \code{\linkS4class{CHMfactor}}
  (rather than from class \code{Cholesky}) that specify the Cholesky
  factor in \emph{sparse} (rather than dense) format.
}
\seealso{
  Classes \code{\linkS4class{dtrMatrix}} and
  \code{\linkS4class{dpoMatrix}}.
  
  Generic functions \code{\link{chol}} and \code{\link{Cholesky}},
  but see \sQuote{Note}.
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
\keyword{classes}
\keyword{algebra}
