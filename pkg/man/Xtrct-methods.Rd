\name{[-methods}
\docType{methods}
\alias{[-methods}
\title{Methods for "[": Extraction or Subsetting in Package 'Matrix'}
%
\alias{[,Matrix,missing,missing,ANY-method}
\alias{[,Matrix,missing,vector,missing-method}
\alias{[,Matrix,vector,missing,missing-method}
\alias{[,Matrix,vector,vector,missing-method}
\alias{[,Matrix,ANY,ANY,ANY-method}% bail out if no explicit method
%
\alias{[,denseMatrix,vector,missing,logical-method}
\alias{[,denseMatrix,vector,vector,logical-method}
\alias{[,denseMatrix,missing,vector,logical-method}
%
\alias{[,gTMatrix,numeric,missing,logical-method}
\alias{[,gTMatrix,missing,numeric,logical-method}
\alias{[,gTMatrix,numeric,numeric,logical-method}
%
\alias{[,dgTMatrix,missing,missing,ANY-method}
\alias{[,dgTMatrix,missing,numeric,logical-method}
\alias{[,dgTMatrix,missing,numeric,missing-method}
\alias{[,dgTMatrix,numeric,missing,logical-method}
\alias{[,dgTMatrix,numeric,missing,missing-method}
\alias{[,dgTMatrix,numeric,numeric,logical-method}
\alias{[,dgTMatrix,numeric,numeric,missing-method}
%
\alias{[,dsparseMatrix,missing,numeric,logical-method}
\alias{[,dsparseMatrix,numeric,missing,logical-method}
\alias{[,dsparseMatrix,numeric,numeric,logical-method}
%
\alias{[,lsparseMatrix,missing,numeric,logical-method}
\alias{[,lsparseMatrix,numeric,missing,logical-method}
\alias{[,lsparseMatrix,numeric,numeric,logical-method}
%-------
\description{
  Methods for \code{"["}, i.e., extraction or subsetting mostly of
  matrices, in package \pkg{Matrix}.
}
\section{Methods}{
  There are more than these:
  \describe{
    \item{x = "Matrix", i = "missing", j = "missing", drop= "ANY"}{ ... }
    \item{x = "Matrix", i = "numeric", j = "missing", drop= "missing"}{ ... }
    \item{x = "Matrix", i = "missing", j = "numeric", drop= "missing"}{ ... }

    \item{x = "dsparseMatrix", i = "missing", j = "numeric", drop= "logical"}{ ... }
    \item{x = "dsparseMatrix", i = "numeric", j = "missing", drop= "logical"}{ ... }
    \item{x = "dsparseMatrix", i = "numeric", j = "numeric", drop= "logical"}{ ... }
  }
}
\seealso{
  \code{\link{[<--methods}} for sub\emph{assign}ment to \code{"Matrix"}
  objects.
  \code{\link{Extract}} about the standard extraction.
}
\examples{
% regression tests are in ../tests/indexing.R
str(m <- Matrix(round(rnorm(7*4),2), nrow = 7))
stopifnot(identical(m, m[]))
m[2, 3]   # simple number
m[2, 3:4] # simple numeric of length 2
m[2, 3:4, drop=FALSE] # sub matrix of class 'dgeMatrix'
## rows or columns only:
m[1,]     # first row, as simple numeric vector
m[,1:2]   # sub matrix of first two columns
}
\keyword{methods}
\keyword{array}
