\name{isSymmetric-methods}
\title{Test if a Matrix is Symmetric}
%
\docType{methods}
\keyword{array}
\keyword{programming}
\keyword{methods}
%
\alias{isSymmetric}
\alias{isSymmetric-methods}
%
\alias{isSymmetric,CsparseMatrix-method}
\alias{isSymmetric,RsparseMatrix-method}
\alias{isSymmetric,TsparseMatrix-method}
\alias{isSymmetric,denseMatrix-method}
\alias{isSymmetric,diagonalMatrix-method}
\alias{isSymmetric,indMatrix-method}
% tolerating numerical fuzz for [dz]Matrix:
\alias{isSymmetric,dgCMatrix-method}
\alias{isSymmetric,dgRMatrix-method}
\alias{isSymmetric,dgTMatrix-method}
\alias{isSymmetric,dgeMatrix-method}
\alias{isSymmetric,dtCMatrix-method}
\alias{isSymmetric,dtRMatrix-method}
\alias{isSymmetric,dtTMatrix-method}
\alias{isSymmetric,dtpMatrix-method}
\alias{isSymmetric,dtrMatrix-method}
%
\description{
Tests whether its argument is a symmetric square matrix, by default
tolerating some numerical fuzz and requiring symmetric \code{dimnames}
in addition to symmetry in the mathematical sense.

\code{isSymmetric} is a generic function in \pkg{base}, which has a
\link[=isSymmetric]{method} for traditional matrices of implicit class
\code{matrix}.  Documented here are methods in package \pkg{Matrix}
for subclasses of virtual class \code{\linkS4class{Matrix}}.
}
\usage{
\S4method{isSymmetric}{denseMatrix}(object, trans = "C", checkDN = TRUE, \dots)
\S4method{isSymmetric}{CsparseMatrix}(object, trans = "C", checkDN = TRUE, \dots)
\S4method{isSymmetric}{RsparseMatrix}(object, trans = "C", checkDN = TRUE, \dots)
\S4method{isSymmetric}{TsparseMatrix}(object, trans = "C", checkDN = TRUE, \dots)
\S4method{isSymmetric}{diagonalMatrix}(object, trans = "C", checkDN = TRUE, \dots)
\S4method{isSymmetric}{indMatrix}(object, checkDN = TRUE, \dots)
}
\arguments{
\item{object}{a \code{\linkS4class{Matrix}}.}
\item{trans}{a character string, \code{"C"} or \code{"T"},
  indicating either the conjugate transpose or transpose,
  i.e., whether to test for Hermitian or symmetric \code{object}.
  It matters only in the complex case.  The default is \code{"C"},
  for consistency with \code{\link{isSymmetric.matrix}}.}
\item{checkDN}{a logical indicating whether symmetry of the
  \code{Dimnames} slot of \code{object} should be checked.}
\item{\dots}{further arguments passed from or to other methods.
  Notably, for numeric or complex (but not integer, logical,
  or pattern) \code{object}, these are used to set a tolerance
  for asymmetry.  See arguments \code{tol}, \code{tol1}, and
  \code{\dots} of \code{\link{isSymmetric.matrix}}.}
}
\value{
A logical, either \code{TRUE} or \code{FALSE} (never \code{NA}).
}
\details{
The \code{Dimnames} \link{slot} of \code{object}, say \code{dn},
is considered to be symmetric if and only if
\itemize{
  \item \code{dn[[1]]} and \code{dn[[2]]} are identical \emph{or}
  one is \code{NULL}; \emph{and}
  \item \code{ndn <- names(dn)} is \code{NULL} \emph{or}
  \code{ndn[1]} and \code{ndn[2]} are identical \emph{or}
  one is the empty string \code{""}.
}
Hence \code{list(a=nms, a=nms)} is considered to be symmetric
and so too are \code{list(a=nms, NULL)} and \code{list(NULL, a=nms)}.

Note that this definition is \emph{looser} than that employed by
\code{\link{isSymmetric.matrix}}, which requires \code{dn[1]} and
\code{dn[2]} to be identical, where \code{dn} is the \code{dimnames}
attribute of a traditional matrix.
}
\seealso{
\code{\link{forceSymmetric}}, \code{\link{symmpart}}, \code{\link{skewpart}};
virtual class \code{\linkS4class{symmetricMatrix}} and its subclasses.
}
\examples{
isSymmetric(Diagonal(4)) # TRUE of course
M <- Matrix(c(1,2,2,1), 2,2)
isSymmetric(M) # TRUE (*and* of formal class "dsyMatrix")
isSymmetric(as(M, "generalMatrix")) # still symmetric, even if not "formally"
isSymmetric(triu(M)) # FALSE

## Look at implementations:
showMethods("isSymmetric", includeDefs = TRUE) # includes S3 generic from base
}
