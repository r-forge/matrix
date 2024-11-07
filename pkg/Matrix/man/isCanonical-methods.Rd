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
isCanonical(x, \dots)
}
\arguments{
  \item{x}{an \R{} object.}
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
## NA->TRUE for nonzero pattern
m1 <- m1. <- new("ngeMatrix", Dim = c(2L, 2L),
                 x = c(FALSE, NA, TRUE, NA))
m1@x <- c(FALSE, TRUE, TRUE, TRUE)
m1; m1. # same class, same matrix, different representations

## transpose for symmetric, conjugate transpose for Hermitian
m2 <- m2. <- new("dsyMatrix", Dim = c(2L, 2L), uplo = "U",
                 x = c(1, 2, 4, 8))
m2@x <- c(1, 4, 4, 8)
m2; m2.

## zero for triangular, always with diag="U"->"N"
m3 <- m3. <- new("dtrMatrix", Dim = c(2L, 2L), uplo = "L", diag = "U",
                 x = c(1, 2, 4, 8))
m3@x <- c(1, 2, 0, 1); m3@diag <- "N"
m3; m3.

## nothing to zero in sparse cases, but unit diagonal must be allocated
m4 <- m4. <- new("ltCMatrix", Dim = c(2L, 2L), uplo = "U", diag = "U",
                 p = c(0L, 0L, 0L))
m4@p <- 0:2; m4@i <- 0:1; m4@x <- c(TRUE, TRUE); m4@diag <- "N"
m4; m4.

## similarly
m5 <- m5. <- new("ndiMatrix", Dim = c(2L, 2L), diag = "U")
m5@x <- c(TRUE, TRUE); m5@diag <- "N"
m5; m5.

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

nms <- paste0(rep(c("m", "v"), c(5L, 2L)), c(1:5, 1:2))
obj <- mget(nms, inherits = FALSE)
obj. <- mget(paste0(nms, "."), inherits = FALSE)
stopifnot( all(vapply(obj , isCanonical, FALSE)),
          !any(vapply(obj., isCanonical, FALSE)))
}
