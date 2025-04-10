\name{bind-methods}
\title{'cbind()' and 'rbind()' recursively built on cbind2/rbind2}
%
\docType{methods}
\keyword{array}
\keyword{manip}
\keyword{methods}
%
\alias{cbind2}
\alias{cbind2-methods}
\alias{rbind2}
\alias{rbind2-methods}
\alias{bind-methods}
%
\alias{cbind2,Matrix,Matrix-method}
\alias{cbind2,Matrix,NULL-method}
\alias{cbind2,Matrix,matrix-method}
\alias{cbind2,Matrix,missing-method}
\alias{cbind2,Matrix,vector-method}
\alias{cbind2,NULL,Matrix-method}
\alias{cbind2,matrix,Matrix-method}
\alias{cbind2,vector,Matrix-method}
%
\alias{rbind2,Matrix,Matrix-method}
\alias{rbind2,Matrix,NULL-method}
\alias{rbind2,Matrix,matrix-method}
\alias{rbind2,Matrix,missing-method}
\alias{rbind2,Matrix,vector-method}
\alias{rbind2,NULL,Matrix-method}
\alias{rbind2,matrix,Matrix-method}
\alias{rbind2,vector,Matrix-method}
%
\alias{c.Matrix}
\alias{c.sparseVector}
%
\description{
  The base functions \code{\link{cbind}} and \code{\link{rbind}} are
  defined for an arbitrary number of arguments and hence have the first
  formal argument \code{...}.  Now, when S4 objects are found among the arguments,
  base \code{cbind()} and \code{rbind()} internally \dQuote{dispatch}
  \emph{recursively}, calling \code{\link{cbind2}} or \code{\link{rbind2}}
  respectively, where these have methods defined and so should dispatch
  appropriately.

  \code{\link{cbind2}()} and \code{\link{rbind2}()} are from the
  \pkg{methods} package, i.e., standard \R, and have been provided for
  binding together \emph{two} matrices, where in \pkg{Matrix}, we have
  defined methods for these and the \code{'Matrix'} matrices.
}
\usage{
## cbind(..., deparse.level = 1)
## rbind(..., deparse.level = 1)

\S4method{cbind2}{Matrix,Matrix}(x, y, \dots)
\S4method{rbind2}{Matrix,Matrix}(x, y, \dots)
}
\arguments{
  \item{\dots}{for \code{[cr]bind}, vector- or matrix-like \R objects
	to be bound together; for \code{[cr]bind2}, further arguments
	passed to or from methods; see \code{\link{cbind}} and
	\code{\link[methods]{cbind2}}.}
  \item{deparse.level}{integer controlling the construction of labels
	in the case of non-matrix-like arguments; see \code{\link{cbind}}.}
  \item{x, y}{vector- or matrix-like \R objects to be bound together.}
}
\value{
  typically a \sQuote{matrix-like} object of a similar
  \code{\link{class}} as the first argument in \code{\dots}.

  Note that sometimes by default, the result is a
  \code{\linkS4class{sparseMatrix}} if one of the arguments is (even in
  the case where this is not efficient).  In other cases,
  the result is chosen to be sparse when there are more zero entries is
  than non-zero ones (as the default \code{sparse} in
  \code{\link{Matrix}()}).
}
\author{Martin Maechler}
\seealso{\code{\link{cbind}}, \code{\link[methods]{cbind2}}.

  Our class definition help pages mentioning \code{cbind2()} and
  \code{rbind2()} methods:
  \code{"\linkS4class{denseMatrix}"},
  \code{"\linkS4class{diagonalMatrix}"},
  \code{"\linkS4class{indMatrix}"}.
}
\examples{
(a <- matrix(c(2:1,1:2), 2,2))

(M1 <- cbind(0, rbind(a, 7))) # a traditional matrix

D <- Diagonal(2)
(M2 <- cbind(4, a, D, -1, D, 0)) # a sparse Matrix

stopifnot(validObject(M2), inherits(M2, "sparseMatrix"),
          dim(M2) == c(2,9))
}

%% \item{cbind2}{\code{(x = *, y = *)}: several methods for binding
%%   matrices together, column-wise, see the basic \code{\link{cbind}}
%%   and \code{\link{rbind}} functions.\cr
%%   Note that the result will typically be sparse, even when one
%%   argument is dense and larger than the sparse one.}
%% \item{rbind2}{\code{(x = *, y = *)}: binding matrices together
%%   row-wise, see \code{cbind2} above.}
