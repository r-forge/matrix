\name{isCanonical-methods}
\title{Test if an Object is Canonical}
%
\docType{methods}
\keyword{array}
\keyword{programming}
\keyword{methods}
%
\alias{isCanonical}
\alias{isCanonical-methods}
%
\alias{isCanonical,CsparseMatrix-method}
\alias{isCanonical,RsparseMatrix-method}
\alias{isCanonical,TsparseMatrix-method}
\alias{isCanonical,denseMatrix-method}
\alias{isCanonical,diagonalMatrix-method}
\alias{isCanonical,indMatrix-method}
\alias{isCanonical,sparseVector-method}
%
\description{
  Where the representation of an abstract object by an \R{} object
  is not unique (because some details are unspecified by the class
  and unreferenced by methods), test for a \dQuote{canonical}
  representation.
}
\usage{
isCanonical(object, \dots)
}
\arguments{
  \item{object}{an \R{} object.}
  \item{\dots}{optional arguments passed to methods.}
}
\details{
  What is meant by \dQuote{canonical} is up to the method writer.  In
  particular, a canonical representation need not be unique.  In other
  words, \code{canonicalize(x)} and \code{canonicalize(y)} need not be
  identical for all sufficiently similar \code{x} and \code{y} of a
  given class.
}
%% TODO: document what "canonical" means for our classes
\value{
  \code{TRUE} if \code{x} is canonical and \code{FALSE} otherwise.
}
\seealso{\code{\link{forceCanonical}}.}
\examples{
## NA->TRUE for nonzero pattern dense ("ndenseMatrix")
m1 <- m1. <- new("ngeMatrix", Dim = c(2L, 2L),
                 x = c(FALSE, NA, TRUE, NA))
m1@x <- c(FALSE, TRUE, TRUE, TRUE)
m1; m1. # same class, same matrix, different representations

## transpose for symmetric, conjugate transpose for Hermitian
m2 <- m2. <- new("zsyMatrix", Dim = c(2L, 2L), uplo = "U", trans = "C",
                 x = c(1+1i, 2+2i, 4+4i, 8+8i))
m2@x <- c(1+0i, 4-4i, 4+4i, 8+0i)
m2; m2.

## zero for triangular, setting diagonal if diag="U"
m3 <- m3. <- new("dtrMatrix", Dim = c(2L, 2L), uplo = "L", diag = "U",
                 x = c(1, 2, 4, 8))
m3@x <- c(1, 2, 0, 1)
m3; m3.

## integer-valued (and where possible integer-typed) slots
## for sparse vectors
v1 <- v1. <- new("nsparseVector", length = 4.2, i = c(2.1, 4.5))
v1@length <- 4L; v1@i <- c(2L, 4L)
v1; v1.

## similarly
im <- .Machine[["integer.max"]]
v2 <- v2. <- new("nsparseVector", length = im + 1.5, i = im + 0.5)
v2@length <- im + 1; v2@i <- im
## v2; v2. # don't print such long vectors

nms <- paste0(rep(c("m", "v"), c(5L, 2L)), c(1:3, 1:2))
obj <- mget(nms, inherits = FALSE)
obj. <- mget(paste0(nms, "."), inherits = FALSE)
stopifnot( all(vapply(obj , isCanonical, NA)),
          !any(vapply(obj., isCanonical, NA)))
}
