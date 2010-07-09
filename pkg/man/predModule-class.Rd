\name{predModule-class}
\Rdversion{1.1}
\title{Class "predModule" and SubClasses}
\docType{class}
\alias{predModule-class}
\alias{dPredModule-class}
\alias{sPredModule-class}
%
\alias{coerce,ddenseModelMatrix,predModule-method}
\alias{coerce,dsparseModelMatrix,predModule-method}
\description{
  The class \code{"predModule"} and notably its subclasses
  \code{"dPredModule"} and \code{"sPredModule"} encapsulate information
  about linear predictors in statistical models.  They incorporate a
  \code{\linkS4class{modelMatrix}}, the corresponding coefficients and a
  representation of a triangular factor from the, possibly weighted or
  otherwise modified, model matrix.
}
\section{Objects from the Classes}{
  Objects are typically created by coercion from objects of class
  \code{\linkS4class{ddenseModelMatrix}} or
  \code{\linkS4class{dsparseModelMatrix}}.
}
\section{Slots}{
  The virtual class \code{"predModule"} contains a
  \code{\linkS4class{modelMatrix}} and one additional slot,
  \describe{
    \item{\code{coef}:}{\code{"numeric"} coefficient vector of length
      \code{ncol(.)}.
    }
  }
  The actual classes \code{"dPredModule"} and \code{"sPredModule"}
  contain an additional slot
  \describe{
    \item{\code{fac}:}{a representation of a triangular factor.  For the
      \code{"dpredModule"} class this factor is a
      \code{\linkS4class{Cholesky}} object.  For the \code{"spredModule"} class
      it is a \code{\linkS4class{CHMfactor}} object.
    }
  }
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "ddenseModelMatrix", to = "predModule")}:
      Creates a \code{"dPredModule"} object.}
    \item{coerce}{\code{signature(from = "dsparseModelMatrix", to = "predModule")}:
      Creates an \code{"sPredModule"} object.}
  }
}
\author{Douglas Bates}
\seealso{
  \code{\link{model.Matrix}} which is a simple wrapper around the traditional
  \code{\link{model.matrix}} and returns a \code{"ddenseModelMatrix"}
  object or, optionally, a \code{"dsparseModelMatrix"} object which can
  then be coerced to a sparse or dense \code{"predModule"}.
}
\examples{
showClass("dPredModule")
showClass("sPredModule")

## see   example(model.Matrix)
}
\keyword{classes}
