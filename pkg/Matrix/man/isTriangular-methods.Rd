\name{isTriangular-methods}
\title{Test if a Matrix is Triangular}
%
\docType{methods}
\keyword{array}
\keyword{programming}
\keyword{methods}
%
\alias{isTriangular}
\alias{isTriangular-methods}
%
\alias{isTriangular,CsparseMatrix-method}
\alias{isTriangular,RsparseMatrix-method}
\alias{isTriangular,TsparseMatrix-method}
\alias{isTriangular,denseMatrix-method}
\alias{isTriangular,diagonalMatrix-method}
\alias{isTriangular,indMatrix-method}
\alias{isTriangular,matrix-method}
%
\description{
  \code{isTriangular} tests whether its argument
  is a triangular matrix.  Unlike the analogous
  \code{\link{isSymmetric}}, it is generically
  from \pkg{Matrix} rather than \pkg{base}.  Hence \pkg{Matrix}
  defines methods for traditional matrices of implicit \link{class}
  \code{"\link{matrix}"} in addition to matrices inheriting from
  virtual class \code{"\linkS4class{Matrix}"}.

  By our definition, triangular matrices are \emph{square},
  i.e., they have the same number of rows and columns.
}
\usage{
isTriangular(object, upper = NA, \dots)
}
\arguments{
  \item{object}{an \R object, typically a matrix.}
  \item{upper}{a \link{logical}, either \code{TRUE} or \code{FALSE},
    in which case \code{TRUE} is returned only for upper or lower
    triangular \code{object}; or otherwise \code{\link{NA}} (the
    default), in which case \code{TRUE} is returned for any triangular
    \code{object}.}
  \item{\dots}{optional arguments passed from or to other methods.}
}
\value{
  A \link{logical}, either \code{TRUE} or \code{FALSE}
  (never \code{\link{NA}}).

  If \code{object} is triangular and \code{upper} is \code{NA}, then
  \code{isTriangular} returns \code{TRUE} with an \link[=attr]{attribute}
  \code{kind}, either \code{"U"} or \code{"L"}, indicating that
  \code{object} is \bold{u}pper or \bold{l}ower triangular, respectively.
  Users should not rely on how \code{kind} is determined for diagonal
  matrices, which are both upper and lower triangular.
}
\seealso{
  \code{\link{isSymmetric}};
  virtual class \code{"\linkS4class{triangularMatrix}"} and subclasses.
}
\examples{
isTriangular(Diagonal(4))
## is TRUE: a diagonal matrix is also (both upper and lower) triangular
(M <- Matrix(c(1,2,0,1), 2,2))
isTriangular(M) # TRUE (*and* of formal class "dtrMatrix")
isTriangular(as(M, "generalMatrix")) # still triangular, even if not "formally"
isTriangular(crossprod(M)) # FALSE

## Look at implementations:
showMethods("isTriangular", includeDefs = TRUE)
}
