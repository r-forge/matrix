\name{t-methods}
%%
%% \alias{t}
\alias{t-methods}
\alias{t,CsparseMatrix-method}
\alias{t,RsparseMatrix-method}
\alias{t,TsparseMatrix-method}
\alias{t,denseMatrix-method}
\alias{t,diagonalMatrix-method}
\alias{t,indMatrix-method}
\alias{t,pMatrix-method}
\alias{t,sparseVector-method}
\alias{.tCRT}
%%
\title{Matrix Transpose}
\description{
  S4 methods defined in package \pkg{Matrix} for S4 generic function
  \code{t}, generalizing \code{base::\link{t}}.
}
\usage{
\S4method{t}{CsparseMatrix}(x)
\S4method{t}{RsparseMatrix}(x)
\S4method{t}{TsparseMatrix}(x)
\S4method{t}{denseMatrix}(x)
\S4method{t}{diagonalMatrix}(x)
\S4method{t}{indMatrix}(x)
\S4method{t}{pMatrix}(x)
\S4method{t}{sparseVector}(x)

.tCRT(x, lazy = TRUE)
}
\arguments{
  \item{x}{a \code{\linkS4class{Matrix}} object or
    \code{\linkS4class{sparseVector}} object.  Vectors must have length
    in the range of type \code{\link{integer}}; they are handled as
    one-column matrices.\cr
	\cr
    For \code{.tCRT}, \code{x} must inherit from one of
    \code{[CRT]sparseMatrix}.}
  \item{lazy}{a logical indicating if transposition can change the
    storage format.}
}
\details{
  For consistency, methods for \code{t} preserve the storage format of
  matrix arguments, even if that is not most efficient.  In particular,
  if \code{x} is a \code{\linkS4class{CsparseMatrix}} object,
  then \code{t(x)} is a \code{CsparseMatrix} object, even though the
  \code{\linkS4class{RsparseMatrix}} representation of the transpose
  is obtained more cheaply.  For \dQuote{lazy} transposing of compressed
  sparse column or row format matrices, interchanging the column and row
  formats, use \code{.tCRT}.

  If \code{x} is formally symmetric
  (see \code{\linkS4class{symmetricMatrix}}),
  then \code{x} and \code{t(x)} represent the same matrix with opposite
  \code{uplo} slots, such that one stores only the upper triangle and
  the other stores only the lower triangle.
}
\value{
  A \code{\linkS4class{Matrix}} object representing the transpose of
  \code{x}.  The class is a superclass of the class of \code{x} unless
  \code{x} is a \code{\linkS4class{sparseVector}} object.
}
\seealso{
  \code{\linkS4methods{ct}} and \code{\link{.ctCRT}} for matrix
  conjugate transpose.  \code{\linkS4methods{crossprod}} for matrix
  products involving transposes.
}
\examples{
(X <- new("igRMatrix",
          Dim = c(3L, 5L),
          Dimnames = list(A = paste0("a", 1:3), B = paste0("b", 1:5)),
          p = c(0L, 1L, 2L, 5L), j = c(1L, 0L, 0L, 1L, 4L), x = -2:2))
( tX <-  t   (X)) # compressed sparse row
(.tX <- .tCRT(X)) # compressed sparse column

stopifnot(identical(dim(tX), dim(X)[2:1]),
          identical(dimnames(tX), dimnames(X)[2:1]),
          identical(X,  t   ( tX)),
          identical(X, .tCRT(.tX)),
          identical( tX, as(.tX, "RsparseMatrix")),
          identical(.tX, as( tX, "CsparseMatrix")),
          identical(tX, .tCRT(X, lazy = FALSE)))

(S <- crossprod(X)) # formally symmetric
(tS <- t(S)) # *output* is identical but stored triangle is different
stopifnot(is(S, "dpRMatrix"),
          all.equal(S, tS, tolerance = 0),
          S@uplo == "U", tS@uplo == "L",
          identical(S, t(tS)))
}
\docType{methods}
\keyword{algebra}
\keyword{array}
\keyword{methods}
