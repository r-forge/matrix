\name{BunchKaufman-methods}
\title{Methods for Bunch-Kaufman Factorization}
%
\docType{methods}
\keyword{algebra}
\keyword{array}
\keyword{methods}
%
\alias{BunchKaufman}
\alias{BunchKaufman-methods}
%
\alias{BunchKaufman,dspMatrix-method}
\alias{BunchKaufman,dsyMatrix-method}
\alias{BunchKaufman,matrix-method}
%
\description{
  Computes the Bunch-Kaufman factorization of an \eqn{n \times n}{n-by-n}
  real, symmetric matrix \eqn{A}, which has the general form
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
  
  Methods are built on LAPACK routines \code{DSYTRF} and \code{DSPTRF}.
}
\usage{
BunchKaufman(x, \dots)
\S4method{BunchKaufman}{dsyMatrix}(x, \dots)
\S4method{BunchKaufman}{dspMatrix}(x, \dots)
\S4method{BunchKaufman}{matrix}(x, uplo = "U", \dots)
}
\arguments{
  \item{x}{a \link[=is.finite]{finite} symmetric matrix or
    \code{\linkS4class{Matrix}} to be factorized.
    If \code{x} is square but not symmetric, then it will be
    \emph{treated} as symmetric; see \code{uplo}.}
  \item{uplo}{a string, either \code{"U"} or \code{"L"},
    indicating which triangle of \code{x} should be used
    to compute the factorization.}
  \item{\dots}{further arguments passed to or from methods.}
}
\value{
  An object representing the factorization, inheriting from
  virtual class \code{\linkS4class{BunchKaufmanFactorization}}.
  The specific class is \code{\linkS4class{BunchKaufman}} unless
  \code{x} inherits from virtual class \code{\linkS4class{packedMatrix}},
  in which case it is \code{\linkS4class{pBunchKaufman}}.
}
\seealso{
  Classes \code{\linkS4class{BunchKaufman}} and
  \code{\linkS4class{pBunchKaufman}} and their methods.
  
  Classes \code{\linkS4class{dsyMatrix}} and
  \code{\linkS4class{dspMatrix}}.
  
  Generic function \code{\link{expand2}},
  for constructing matrix factors from the result.
  
  Generic functions \code{\link{Cholesky}}, \code{\link{Schur}},
  \code{\link{lu}}, and \code{\link{qr}},
  for computing other factorizations.
}
\references{
  The LAPACK source code, including documentation; see
  \url{https://netlib.org/lapack/double/dsytrf.f} and
  \url{https://netlib.org/lapack/double/dsptrf.f}.
}
\examples{
showMethods("BunchKaufman", inherited = FALSE)

data(CAex)
dim(CAex)
isSymmetric(CAex) # TRUE

S <- as(CAex, "symmetricMatrix") # dsCMatrix

## No methods for sparseMatrix yet:
if(FALSE)
   BunchKaufman(S)
## But two for *dense*Matrix:
 bk <- BunchKaufman(as(S, "unpackedMatrix"))
pbk <- BunchKaufman(as(S,   "packedMatrix"))
pbk

## See the underlying LAPACK representation
 bk.M <- as( bk, "dtrMatrix")
pbk.M <- as(pbk, "dtpMatrix")
pbk.M[1:6, 1:6]

## See how sparse it is, too
bk.M.s <- as(bk.M, "sparseMatrix")
image(bk.M.s)

## An upper triangular 72-by-72 matrix with only 144 nonzero entries
str(bk.M.s)

stopifnot(is(bk.M,   "triangularMatrix"),
          is(bk.M,     "unpackedMatrix"),
          is(bk.M.s, "triangularMatrix"),
          is(bk.M.s,    "CsparseMatrix"))
}
