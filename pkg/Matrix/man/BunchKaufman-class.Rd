\name{BunchKaufman-class}
\title{Dense Bunch-Kaufman Factorizations}
%
\docType{class}
\keyword{algebra}
\keyword{array}
\keyword{classes}
%
\alias{BunchKaufman-class}
\alias{pBunchKaufman-class}
%
\alias{coerce,BunchKaufman,Matrix-method}
\alias{coerce,pBunchKaufman,Matrix-method}
%
\description{
  Classes \code{BunchKaufman} and \code{pBunchKaufman} represent
  Bunch-Kaufman factorizations of real, symmetric matrices.  They
  store the nonzero entries of multiple \code{n}-by-\code{n} sparse
  matrix factors in dense format as a vector of length \code{n*n}
  (\code{BunchKaufman}) or \code{n*(n+1)/2} (\code{pBunchKaufman}),
  the latter giving the \dQuote{\bold{p}acked} representation.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}}{inherited from virtual class
      \code{\linkS4class{MatrixFactorization}}.}
    \item{\code{uplo}}{a string, either \code{"U"} or \code{"L"},
      indicating which triangle (upper or lower) of the factorized
      symmetric matrix was used to compute the factorization and
      in turn how the \code{x} slot is partitioned.}
    \item{\code{x}}{a numeric vector of length
      \code{n*n} (\code{BunchKaufman}) or
      \code{n*(n+1)/2} (\code{pBunchKaufman}), where \code{n=Dim[1]}.
      The details of the representation are specified by the LAPACK
      manual for its routines \code{DSYTRF} and \code{DSPTRF}.}
    \item{\code{perm}}{an integer vector of length \code{n=Dim[1]}
      specifying row and column interchanges and in turn how the
      \code{x} slot is partitioned.}
  }
}
\section{Extends}{
  Class \code{\linkS4class{BunchKaufmanFactorization}}, directly.
  Class \code{\linkS4class{MatrixFactorization}}, by class
  \code{\linkS4class{BunchKaufmanFactorization}}, distance 2.
}
\section{Instantiation}{
  Objects can be generated directly by calls of the form
  \code{new("BunchKaufman", ...)} or \code{new("pBunchKaufman", ...)},
  but they are more typically obtained as the value of
  \code{\link{BunchKaufman}(x)} for \code{x} inheriting from
  \code{\linkS4class{dsyMatrix}} or \code{\linkS4class{dspMatrix}}.
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "p?BunchKaufman", to = "Matrix")}:
      returns a \code{dtrMatrix} or \code{dtpMatrix}, useful for
      inspecting the internal representation of the factorization;
      see Note.}
  }
}
\note{
  In \pkg{Matrix} \code{< 1.5-5}, class \code{BunchKaufman} extended
  \code{\linkS4class{dtrMatrix}} and class \code{pBunchKaufman}
  extended \code{\linkS4class{dtpMatrix}}, reflecting the fact that
  the internal representation of the factorization is fundamentally
  triangular: there are \code{n*(n+1)/2} \dQuote{parameters}, and
  these can be arranged systematically to form an \code{n}-by-\code{n}
  triangular matrix.

  \pkg{Matrix} 1.5-5 removed these extensions so that
  \code{BunchKaufman} and \code{pBunchKaufman} would no longer
  inherit methods from \code{dtrMatrix} and \code{dtpMatrix}.
  The availability of such methods gave the wrong impression that
  \code{BunchKaufman} and \code{pBunchKaufman} represent proper
  dense triangular matrices, when in fact they represent multiple
  sparse matrix factors, not all triangular, in a dense format.

  The coercions \code{as(., "Matrix")} are provided for users who
  understand these caveats.
}
\seealso{
  Class \code{\linkS4class{dsyMatrix}} and its packed counterpart.

  Generic function \code{\link{BunchKaufman}}.
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
