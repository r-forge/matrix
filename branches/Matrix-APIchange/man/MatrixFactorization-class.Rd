\name{MatrixFactorization-class}
\alias{MatrixFactorization-class}
\docType{class}
\title{Class "MatrixFactorization" of Matrix Factorizations}
%
\alias{expand,MatrixFactorization-method}
\alias{show,MatrixFactorization-method}
\description{
  The class \code{"MatrixFactorization"} is the virtual (super) class of
  (potentially) all matrix factorizations of matrices from package
  \pkg{Matrix}.
}
\section{Objects from the Class}{A virtual Class: No objects may be
  created from it.}
\section{Slots}{
  \describe{
    \item{\code{Dim}:}{Object of class \code{"integer"} - the dimensions
      of the original matrix - must be an integer vector with exactly two
      non-negative values.}
  }
}
\section{Methods}{
  \describe{
    \item{expand}{\code{signature(x = "MatrixFactorization")}: this has
      not been implemented yet for all matrix factorizations. It should
      return a list whose components are matrices which when multiplied
      return the original \code{\linkS4class{Matrix}} object.}
    \item{show}{\code{signature(object = "MatrixFactorization")}: ... }
	 }
}
\seealso{
  classes inheriting from \code{"MatrixFactorization"}, such as
  \code{\linkS4class{CHMfactor}},
  \code{\linkS4class{LU}}, and
  \code{\linkS4class{sparseQR}}.
}
\examples{
showClass("MatrixFactorization")
}
\keyword{classes}
