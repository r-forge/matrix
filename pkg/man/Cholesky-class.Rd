\name{Cholesky-class}
\docType{class}
\alias{Cholesky-class}

\title{Cholesky decompositions}
\description{The \code{"Cholesky"} class is the class of Cholesky
  decomposition of a positive-semidefinite, real matrices}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("Cholesky", ...)}
  or by calls of the form \code{chol(pm)} where \code{pm} inherits from
  the \code{"dpoMatrix"} class or as a side-effect of other functions
  applied to \code{"dpoMatrix"} objects.
}
\section{Slots}{
  \describe{
    \item{\code{uplo}:}{inherited from the \code{"dtrMatrix"} class.}
    \item{\code{diag}:}{inherited from the \code{"dtrMatrix"} class.}
    \item{\code{x}:}{inherited from the \code{"dtrMatrix"} class.}
    \item{\code{Dim}:}{inherited from the \code{"dtrMatrix"} class.}
    \item{\code{rcond}:}{inherited from the \code{"dtrMatrix"} class.}
    \item{\code{factors}:}{inherited from the \code{"dtrMatrix"}
      class.  For this class the \code{factors} slot is always an
      empty list.}
  }
}
\section{Extends}{
Class \code{"dtrMatrix"}, directly.
Class \code{"dgeMatrix"}, by class \code{"dtrMatrix"}.
Class \code{"Matrix"}, by class \code{"dtrMatrix"}.
}
\section{Methods}{
No methods defined with class "Cholesky" in the signature.
}
\seealso{
  \code{\link{dtrMatrix-class}}, \code{chol}, \code{dpoMatrix-class}
}
\examples{
print(pm <- crossprod(Matrix(rnorm(18), nrow = 6, ncol = 3)))
print(ch <- chol(pm))
if (toupper(ch@uplo) == "U") crossprod(ch)
}
\keyword{classes}
\keyword{algebra}
