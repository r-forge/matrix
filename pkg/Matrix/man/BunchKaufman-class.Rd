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
\alias{coerce,BunchKaufman,dtrMatrix-method}
\alias{coerce,pBunchKaufman,dtpMatrix-method}
\alias{determinant,BunchKaufman,logical-method}
\alias{determinant,pBunchKaufman,logical-method}
%
\description{
  Classes \code{BunchKaufman} and \code{pBunchKaufman} represent
  Bunch-Kaufman factorizations of \eqn{n \times n}{n-by-n} real,
  symmetric matrices \eqn{A}, having the general form
  \deqn{A = L D_{L} L' = U' D_{U} U}{A = L * DL * L' = U' * DU * U}
  where
  \eqn{D_{L}}{DL} and \eqn{D_{U}}{DU} are symmetric, block diagonal
  matrices composed of \eqn{b_{L}}{bL} and \eqn{b_{U}}{bU}
  \eqn{1 \times 1}{1-by-1} or \eqn{2 \times 2}{2-by-2} diagonal blocks,
  and
  \eqn{L = \prod_{k = 1}^{b_{L}} P_{k} L_{k}}{prod(Pk * Lk : k = 1,...,bL)}
  is the product of \eqn{b_{L}}{bL} row-permuted unit lower triangular
  matrices, each having nonzero entries below the diagonal in 1 or 2 columns,
  and 
  \eqn{U = \prod_{k = 1}^{b_{U}} P_{k} U_{k}}{prod(Pk * Uk : k = 1,...,bU)}
  is the product of \eqn{b_{U}}{bU} row-permuted unit upper triangular
  matrices, each having nonzero entries above the diagonal in 1 or 2 columns.
  
  These classes store the nonzero entries of the
  \eqn{2 b_{L} + 1}{2*bL+1} or \eqn{2 b_{U} + 1}{2*bU+1} factors,
  which are individually sparse,
  in a dense format as a vector of length
  \eqn{nn}{n*n} (\code{BunchKaufman}) or
  \eqn{n(n+1)/2}{n*(n+1)/2} (\code{pBunchKaufman}),
  the latter giving the \dQuote{packed} representation.
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
      The details of the representation are specified by the manual
      for LAPACK routines \code{DSYTRF} and \code{DSPTRF}.}
    \item{\code{perm}}{an integer vector of length \code{n=Dim[1]}
      specifying row and column interchanges as described in the manual
      for LAPACK routines \code{DSYTRF} and \code{DSPTRF}.}
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
    \item{\code{coerce}}{\code{signature(from = "BunchKaufman", to = "dtrMatrix")}:
      returns a \code{\linkS4class{dtrMatrix}}, useful for inspecting
      the internal representation of the factorization; see \sQuote{Note}.}
    \item{\code{coerce}}{\code{signature(from = "pBunchKaufman", to = "dtpMatrix")}:
      returns a \code{\linkS4class{dtpMatrix}}, useful for inspecting
      the internal representation of the factorization; see \sQuote{Note}.}
    \item{\code{determinant}}{\code{signature(from = "p?BunchKaufman", logarithm = "logical")}:
      computes the determinant of the factorized matrix \eqn{A}
      or its logarithm.}
    \item{\code{expand2}}{\code{signature(x = "p?BunchKaufman")}:
      see \code{\link{expand2-methods}}.}
  }
}
\note{
  In \pkg{Matrix} \code{< 1.5-5}, class \code{BunchKaufman} extended
  \code{\linkS4class{dtrMatrix}} and class \code{pBunchKaufman} extended
  \code{\linkS4class{dtpMatrix}}, reflecting the fact that the internal
  representation of the factorization is fundamentally triangular:
  there are \code{n*(n+1)/2} \dQuote{parameters}, and these can be
  arranged systematically to form a \code{n}-by-\code{n} triangular
  matrix.

  \pkg{Matrix} 1.5-5 removed these extensions so that
  \code{BunchKaufman} and \code{pBunchKaufman} would no longer inherit
  methods from \code{dtrMatrix} and \code{dtpMatrix}.  The availability
  of such methods gave the wrong impression that \code{BunchKaufman}
  and \code{pBunchKaufman} represent proper dense triangular matrices,
  when in fact they represent multiple sparse matrix factors, not all
  triangular, in a dense format.
  
  The coercions \code{as(., "triangularMatrix")} are provided for users
  who understand these caveats.
}
\seealso{
  Class \code{\linkS4class{dsyMatrix}} and its packed counterpart.

  Generic function \code{\link{BunchKaufman}}.
}
\references{
  The LAPACK source code, including documentation; see
  \url{https://netlib.org/lapack/double/dsytrf.f} and
  \url{https://netlib.org/lapack/double/dsptrf.f}.
}
\examples{
(sy <- crossprod(Matrix(rnorm(18L), 6L, 3L)))
(bksy <- BunchKaufman(sy))

(sp <- pack(Matrix(diag(5L) + 1)))
(bksp <- BunchKaufman(sp))

## See the underlying LAPACK representation 
as(bksy, "dtrMatrix")
as(bksp, "dtpMatrix")
}
