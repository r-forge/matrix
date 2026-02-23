\name{which-methods}
%%
%% \alias{which}
\alias{which-methods}
\alias{which,indMatrix-method}
\alias{which,ldenseMatrix-method}
\alias{which,ldiMatrix-method}
\alias{which,lsparseMatrix-method}
\alias{which,lsparseVector-method}
\alias{which,ndenseMatrix-method}
\alias{which,ndiMatrix-method}
\alias{which,nsparseMatrix-method}
\alias{which,nsparseVector-method}
%%
\title{Which Indices are TRUE?}
\description{
  S4 methods defined in package \pkg{Matrix} for S4 generic function
  \code{which}, generalizing \code{base::\link{which}}.
}
\usage{
\S4method{which}{ndenseMatrix}(x, arr.ind = FALSE, useNames = TRUE)
\S4method{which}{ndiMatrix}(x, arr.ind = FALSE, useNames = TRUE)
\S4method{which}{nsparseMatrix}(x, arr.ind = FALSE, useNames = TRUE)
\S4method{which}{nsparseVector}(x, arr.ind = FALSE, useNames = TRUE)

\S4method{which}{indMatrix}(x, arr.ind = FALSE, useNames = TRUE)

\S4method{which}{ldenseMatrix}(x, arr.ind = FALSE, useNames = TRUE)
\S4method{which}{ldiMatrix}(x, arr.ind = FALSE, useNames = TRUE)
\S4method{which}{lsparseMatrix}(x, arr.ind = FALSE, useNames = TRUE)
\S4method{which}{lsparseVector}(x, arr.ind = FALSE, useNames = TRUE)
}
\arguments{
  \item{x}{a \code{\linkS4class{Matrix}} object or
    \code{\linkS4class{sparseVector}} object with pattern or logical
    data.}
  \item{arr.ind}{a logical indicating if array indices are computed and
    returned when \code{x} is a \code{\linkS4class{Matrix}} object.}
  \item{useNames}{a logical, passed to \code{\link{arrayInd}} when
    \code{x} is a \code{\linkS4class{Matrix}} and \code{arr.ind} is
    \code{TRUE}.}
}
\value{
  An integer or double vector or a two-column integer matrix, following
  \code{base::\link{which}}.
}
\examples{
## Column index matrix: matrix whose columns are standard unit vectors
n <- 6L
r <- 4L
(J <- new("indMatrix", Dim = c(n, r), Dimnames = list(u=NULL, v=NULL),
          margin = 2L, perm = sample(n, r, replace=TRUE)))
(w0 <- which(J))
(w1 <- which(J, arr.ind = TRUE))
stopifnot(identical(w0, 0L:(r - 1L) * n + J@perm),
          identical(w1, cbind(u=J@perm, v=1L:r)))
}
\docType{methods}
\keyword{logic}
\keyword{methods}
