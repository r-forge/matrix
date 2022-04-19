\name{isSymmetric-methods}
\title{Methods for Function isSymmetric in Package 'Matrix'}
\docType{methods}
%
\alias{isSymmetric-methods}
\alias{isSymmetric,symmetricMatrix-method}
\alias{isSymmetric,triangularMatrix-method}
\alias{isSymmetric,unpackedMatrix-method}
\alias{isSymmetric,dgeMatrix-method}
\alias{isSymmetric,dtrMatrix-method}
\alias{isSymmetric,lgeMatrix-method}
\alias{isSymmetric,ltrMatrix-method}
\alias{isSymmetric,ngeMatrix-method}
\alias{isSymmetric,ntrMatrix-method}
\alias{isSymmetric,packedMatrix-method}
\alias{isSymmetric,dtpMatrix-method}
\alias{isSymmetric,ltpMatrix-method}
\alias{isSymmetric,ntpMatrix-method}
\alias{isSymmetric,sparseMatrix-method}
\alias{isSymmetric,diagonalMatrix-method}
\alias{isSymmetric,matrix-method}
%
\description{
  \code{isSymmetric(M)} returns a \code{\link{logical}} indicating
  if \code{M} is a symmetric matrix.  This (now) is a \pkg{base} function
  with a default method for the traditional matrices of \code{\link{class}}
  \code{"matrix"}.  Methods here are defined for virtual Matrix classes
  such that it works for all objects inheriting from class
  \code{\linkS4class{Matrix}}.
}
\seealso{\code{\link{forceSymmetric}}, \code{\link{symmpart}}, and the
  formal class (and subclasses) \code{"\linkS4class{symmetricMatrix}"}.
}
%% \usage{
%% ## some methods have a 'tol' :
%% \S4method{isSymmetric}{denseMatrix}(object, tol = 100 * .Machine$double.eps, ...) 
%% \S4method{isSymmetric}{sparseMatrix}(object, tol = 100 * .Machine$double.eps, ...) 
%% }
%% \section{Methods}{
%%   \describe{
%%     \item{object = "symmetricMatrix"}{ ... }
%%     \item{object = "triangularMatrix"}{ ... }
%%     \item{object = "denseMatrix"}{ ... }
%%     \item{object = "diagonalMatrix"}{ ... }
%%     \item{object = "sparseMatrix"}{ ... }
%%   }
%% }
\examples{
isSymmetric(Diagonal(4)) # TRUE of course
M <- Matrix(c(1,2,2,1), 2,2)
isSymmetric(M) # TRUE (*and* of formal class "dsyMatrix")
isSymmetric(as(M, "dgeMatrix")) # still symmetric, even if not "formally"
isSymmetric(triu(M)) # FALSE

## Look at implementations:
showMethods("isSymmetric", includeDefs=TRUE)# "ANY": base's S3 generic; 6 more
}
\keyword{methods}
