\name{forceDiagonal-methods}
\title{Create \code{diagonalMatrix} From the Diagonal Part of a Matrix}
%
\docType{methods}
\keyword{array}
\keyword{methods}
%
\alias{forceDiagonal}
\alias{forceDiagonal-methods}
%
\alias{forceDiagonal,CsparseMatrix-method}
\alias{forceDiagonal,RsparseMatrix-method}
\alias{forceDiagonal,TsparseMatrix-method}
\alias{forceDiagonal,denseMatrix-method}
\alias{forceDiagonal,diagonalMatrix-method}
\alias{forceDiagonal,indMatrix-method}
\alias{forceDiagonal,matrix-method}
%
\description{
   From a \code{matrix} or \code{\link{Matrix}} \code{x}, use the diagonal part,
   \code{\link{diag}(x)} to create a corresponding \code{\linkS4class{diagonalMatrix}}.
}
\usage{
forceDiagonal(x, \dots)

\S4method{forceDiagonal}{matrix}(x, diag = NULL, \dots)
}
\arguments{
  \item{x}{a \code{matrix} or \code{\linkS4class{Matrix}}.}
  \item{diag}{a string, \code{"U"} or \code{"N"} or \code{\link{NULL}},
    determining the \code{diag} slot of the result; if \code{NULL}, as by
    default, the code checks for unit diagonal and uses \code{"U"} in that
    case.}
  \item{\dots}{optional arguments passed from or to other methods.}
}
\value{a \code{\link{Matrix}} inheriting from
  \code{"\linkS4class{diagonalMatrix}"}, even when \code{x} is not diagonal.
}
\seealso{\code{as(x, "diagonalMatrix")} does \emph{not} force, but signal
  an error when \code{x} is not diagonal, see \code{\link{isDiagonal}}.
}
\examples{
forceDiagonal(diag(3)) # --> <ddiMatrix>
forceDiagonal(matrix(1:4, 2))# a <idiMatrix>
stopifnot(identical(forceDiagonal(Matrix(1:4, 2)),
                    Diagonal(x = c(1, 4))))
}
