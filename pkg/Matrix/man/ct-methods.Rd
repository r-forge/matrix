\name{ct-methods}
%%
\alias{ct}
\alias{ct-methods}
\alias{ct,ANY-method}
\alias{ct,CsparseMatrix-method}
\alias{ct,RsparseMatrix-method}
\alias{ct,TsparseMatrix-method}
\alias{ct,denseMatrix-method}
\alias{ct,diagonalMatrix-method}
\alias{ct,indMatrix-method}
\alias{ct,pMatrix-method}
\alias{ct,sparseVector-method}
\alias{.ctCRT}
%%
\title{Matrix Conjugate Transpose}
\description{
  Returns the conjugate transpose of a matrix or data frame, defined as
  \code{\link{Conj}(\link{t}(x))} if \code{x} has complex data and as
  \code{\link{t}(x)} otherwise.
}
\usage{
ct(x)

\S4method{ct}{ANY}(x)
\S4method{ct}{CsparseMatrix}(x)
\S4method{ct}{RsparseMatrix}(x)
\S4method{ct}{TsparseMatrix}(x)
\S4method{ct}{denseMatrix}(x)
\S4method{ct}{diagonalMatrix}(x)
\S4method{ct}{indMatrix}(x)
\S4method{ct}{pMatrix}(x)
\S4method{ct}{sparseVector}(x)

.ctCRT(x, lazy = TRUE)
}
\arguments{
  \item{x}{a \R{} object representing a matrix or vector.
	\code{\linkS4class{Matrix}} object or
    \code{\linkS4class{sparseVector}} object.  Vectors must have length
    in the range of type \code{\link{integer}}; they are handled as
    one-column matrices.\cr
	\cr
    For \code{.ctCRT}, \code{x} must inherit from one of
    \code{[CRT]sparseMatrix}.}
  \item{lazy}{a logical indicating if transposition can change the
    storage format.}
}
\details{
  For consistency, methods for \code{ct} preserve the storage format of
  matrix arguments, even if that is not most efficient.  In particular,
  if \code{x} is a \code{\linkS4class{CsparseMatrix}} object,
  then \code{ct(x)} is a \code{CsparseMatrix} object, even though the
  \code{\linkS4class{RsparseMatrix}} representation of the conjugate
  transpose is obtained more cheaply.  For \dQuote{lazy} transposing of
  compressed sparse column or row format matrices, interchanging the
  column and row formats, use \code{.ctCRT}.

  If \code{x} is formally Hermitian
  (see \code{\linkS4class{symmetricMatrix}}),
  then \code{x} and \code{ct(x)} represent the same matrix with opposite
  \code{uplo} slots, such that one stores only the upper triangle and
  the other stores only the lower triangle.
}
\value{
  A \code{\linkS4class{Matrix}} object representing the conjugate
  transpose of \code{x}.  The class is a superclass of the class of
  \code{x} unless \code{x} is a \code{\linkS4class{sparseVector}}
  object.
}
\seealso{
  \code{\linkS4methods{t}} and \code{\link{.tCRT}} for matrix transpose
  without conjugation.  \code{\linkS4methods{crossprod}} for matrix
  products involving transposes.
}
\examples{
ct(c(1+1i, 2-2i)) # as Conj(t(.)) for complex data

typeof(x2 <- c(1L, 2L))
typeof(print(      t(x2) )) # "integer"
typeof(print(     ct(x2) )) # "integer"
typeof(print(Conj( t(x2)))) # "double" as 'Conj' coerces

stopifnot(identical(ct(1:2), t(1:2)), # keeping non-complex data type
          identical(Conj(t(1:2)), 0 + t(1:2))) 

rcplx <-
function (n)
complex(modulus = rlnorm(n), argument = runif(n, 0, 2 * pi))

set.seed(14L)
(X <- new("zgRMatrix",
          Dim = c(3L, 5L),
          Dimnames = list(A = paste0("a", 1:3), B = paste0("b", 1:5)),
          p = c(0L, 1L, 2L, 5L), j = c(1L, 0L, 0L, 1L, 4L),
          x = round(rcplx(5L), 2L)))
( ctX <-  ct   (X)) # compressed sparse row
(.ctX <- .ctCRT(X)) # compressed sparse column

stopifnot(identical(dim(ctX), dim(X)[2:1]),
          identical(dimnames(ctX), dimnames(X)[2:1]),
          identical(X,  ct   ( ctX)),
          identical(X, .ctCRT(.ctX)),
          identical( ctX, as(.ctX, "RsparseMatrix")),
          identical(.ctX, as( ctX, "CsparseMatrix")),
          identical(ctX, .ctCRT(X, lazy = FALSE)))

(S <- crossprod(X, trans = "C")) # formally Hermitian
(ctS <- ct(S)) # *output* is identical but stored triangle is different
stopifnot(is(S, "zpRMatrix"),
          all.equal(S, ctS, tolerance = 0),
          S@uplo == "U", ctS@uplo == "L",
          identical(S, ct(ctS)))
}
\docType{methods}
\keyword{algebra}
\keyword{array}
\keyword{methods}
