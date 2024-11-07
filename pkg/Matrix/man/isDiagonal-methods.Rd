\name{isDiagonal-methods}
\title{Test if a Matrix is Diagonal}
%
\docType{methods}
\keyword{array}
\keyword{programming}
\keyword{methods}
%
\alias{isDiagonal}
\alias{isDiagonal-methods}
%
\alias{isDiagonal,CsparseMatrix-method}
\alias{isDiagonal,RsparseMatrix-method}
\alias{isDiagonal,TsparseMatrix-method}
\alias{isDiagonal,denseMatrix-method}
\alias{isDiagonal,diagonalMatrix-method}
\alias{isDiagonal,indMatrix-method}
\alias{isDiagonal,matrix-method}
%
\description{
  \code{isDiagonal} tests whether its argument
  is a diagonal matrix.  Unlike the analogous
  \code{\link{isSymmetric}}, it is generically
  from \pkg{Matrix} rather than \pkg{base}.  Hence \pkg{Matrix}
  defines methods for traditional matrices of implicit \link{class}
  \code{"\link{matrix}"} in addition to matrices inheriting from
  virtual class \code{"\linkS4class{Matrix}"}.

  By our definition, diagonal matrices are \emph{square},
  i.e., they have the same number of rows and columns.
}
\usage{
isDiagonal(object, \dots)
}
\arguments{
  \item{object}{an \R object, typically a matrix.}
  \item{\dots}{optional arguments passed from or to other methods.}
}
\value{
  A \link{logical}, either \code{TRUE} or \code{FALSE}
  (never \code{\link{NA}}).
}
\seealso{
  \code{\link{isSymmetric}};
  virtual class \code{"\linkS4class{diagonalMatrix}"} and subclasses.
}
\examples{
isDiagonal(matrix(c(2,0,0,1), 2,2)) # TRUE

## Look at implementations:
showMethods("isDiagonal", includeDefs = TRUE)
}
