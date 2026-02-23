\name{toeplitz-methods}
%%
%% \alias{toeplitz}
\alias{toeplitz-methods}
\alias{toeplitz,sparseVector-method}
%%
\title{Form a Toeplitz Matrix}
\description{
  S4 methods defined in package \pkg{Matrix} for S4 generic function
  \code{toeplitz}, generalizing \code{stats::\link{toeplitz}}.

  The S4 generic function definition is obtained from package
  \pkg{methods} through the \code{\link{implicitGeneric}} mechanism.
}
\usage{
\S4method{toeplitz}{sparseVector}(x, r = NULL, symmetric = is.null(r),
         repr = c("C", "R", "T"), giveCsparse, \dots)
}
\arguments{
  \item{x}{a \code{\linkS4class{sparseVector}} object.}
  \item{r}{a \code{\linkS4class{sparseVector}} object or \code{NULL},
     not used if \code{symmetric} is \code{TRUE}.}
  \item{symmetric}{a logical indicating if the return value is formally
    symmetric.}
  \item{repr}{a character string indicating the format of the return
    value.}
  \item{giveCsparse}{deprecated.}
  \item{\dots}{optional arguments not yet used.}
}
\value{
  A \code{\linkS4class{sparseMatrix}} object representing the Toeplitz
  matrix with first row \code{x} or \code{r} (\code{x} if
  \code{symmetric} is \code{TRUE} or \code{r} is \code{NULL}) and first
  column \code{x}.  \code{x[1]} determines the main diagonal;
  \code{r[1]} is not used.
}
\examples{
x <- new("isparseVector", length = 6L,
         i = c(2L, 4L, 6L), x = -1L:1L)
y <- new("dsparseVector", length = 8L,
         i = c(2L, 8L), x = c(1.1, 2.2))
(t1s <- toeplitz(x))
(t1g <- toeplitz(x, symmetric = FALSE))
(t2g <- toeplitz(x, y))
stopifnot(is(t1s, "isCMatrix"),
          is(t1g, "igCMatrix"),
          is(t2g, "dgCMatrix"),
          identical(t1g, as(t1s, "generalMatrix")),
          identical(tril(t2g[, 1L:6L]), as(tril(t1g), "dMatrix")))
}
\docType{methods}
\keyword{array}
\keyword{ts}
\keyword{methods}
