\name{BunchKaufman-class}
\title{Bunch-Kaufman Factorizations}
\docType{class}
%
\alias{BunchKaufman-class}
\alias{pBunchKaufman-class}
%
\alias{coerce,BunchKaufman,Matrix-method}
\alias{coerce,pBunchKaufman,Matrix-method}
\alias{show,BunchKaufman-method}
\alias{show,pBunchKaufman-method}
%
\description{
  Classes \code{BunchKaufman} and \code{pBunchKaufman} represent
  Bunch-Kaufman factorizations of real, symmetric matrices.  They
  store the nonzero entries of multiple \code{n}-by-\code{n} sparse
  matrix factors in dense format as a vector of length \code{n*n}
  (\code{BunchKaufman}) or \code{n*(n+1)/2} (\code{pBunchKaufman}),
  the latter giving the \dQuote{\bold{p}acked} representation.

  The details of the representation are specified by the LAPACK
  documentation for its routines \code{DSYTRF} and \code{DSPTRF}.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}}{inherited from virtual class
      \code{\linkS4class{MatrixFactorization}}.}
    \item{\code{Dimnames}}{a list of length 2, possibly with
      \link{names}.  Each element must be \code{NULL} or a
      \link{character} vector of length equal to the corresponding
      element of \code{Dim}.}
    \item{\code{uplo}}{a string, either \code{"U"} or \code{"L"},
      indicating which triangle (upper or lower) of the factorized
      symmetric matrix was used to compute the factorization and
      in turn how the \code{x} slot is partitioned.}
    \item{\code{x}}{a vector of type \link{double} and length
      \code{n*n} (\code{BunchKaufman}) or \code{n*(n+1)/2}
      (\code{pBunchKaufman}), if \code{n=Dim[1]}.}
    \item{\code{perm}}{an integer vector of length \code{n=Dim[1]}
      specifying row and column interchanges and in turn how the
      \code{x} slot is partitioned.}
  }
}
\section{Extends}{
  Class \code{\linkS4class{MatrixFactorization}}, directly.
}
\section{Methods}{
  Currently only simple methods for \code{\link{show}} and coercions
  \code{as(x, "Matrix")}.  The latter return \code{dtrMatrix} and
  \code{dtpMatrix} for \code{x} inheriting from \code{Cholesky} and
  \code{pCholesky}, respectively.  The coercions are useful primarily
  for inspecting LAPACK's dense triangular representation of the
  factorization; see \sQuote{Note}.
}
\section{Objects from the class}{
  Objects can be created directly by calls of the form
  \code{new("BunchKaufman", ...)} or \code{new("pBunchKaufman", ...)}.
  More typically, they are obtained as the value of
  \code{\link{BunchKaufman}(x)} for \code{x} inheriting from
  \code{\linkS4class{dsyMatrix}} or \code{\linkS4class{dspMatrix}}.
}
\note{
  In \pkg{Matrix} \code{< 1.5-5}, class \code{BunchKaufman} extended
  \code{\linkS4class{dtrMatrix}} and class \code{pBunchKaufman}
  extended \code{\linkS4class{dtpMatrix}}, reflecting the fact that
  the internal representation of the factorization is fundamentally
  triangular: there are \code{n*(n+1)/2} \dQuote{parameters}, and
  these can be arranged systematically to form an \code{n}-by-\code{n}
  triangular matrix.

  \pkg{Matrix} 1.5-5 removed these \dQuote{extensions} so that
  \code{BunchKaufman} and \code{pBunchKaufman} would \emph{not}
  inherit methods from \code{dtrMatrix} and \code{dtpMatrix}.
  The availability of such methods gave the wrong impression that
  \code{BunchKaufman} and \code{pBunchKaufman} represent proper
  dense triangular matrices, when in fact they represent multiple
  sparse matrix factors, not all triangular, in a dense format.

  The coercions \code{as(., "Matrix")} are provided for users who
  understand these caveats.
}
\seealso{
  Class \code{\linkS4class{dsyMatrix}}.  Generic function
  \code{\link{BunchKaufman}}.
}
\examples{
(sp <- pack(Matrix(diag(5L) + 1)))
(bksp <- BunchKaufman(sp))

(sy <- crossprod(Matrix(rnorm(18L), 6L, 3L)))
(bksy <- BunchKaufman(sy))

## See the underlying LAPACK representation 
as(bksp, "Matrix")
as(bksy, "Matrix")
}
\keyword{classes}
\keyword{algebra}
