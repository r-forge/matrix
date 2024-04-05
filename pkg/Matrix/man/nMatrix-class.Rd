\name{nMatrix-class}
\title{Class "nMatrix" of Non-zero Pattern Matrices}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{nMatrix-class}
%
\alias{coerce,matrix,nMatrix-method}
\alias{coerce,vector,nMatrix-method}
%
\description{
  The \code{nMatrix} class is the virtual \dQuote{mother} class of all
  \emph{\bold{n}on-zero pattern} (or simply \emph{patter\bold{n}})
  matrices in the \pkg{Matrix} package.
}
%\section{Objects from the Class}{A virtual Class: No objects may be
%  created from it.
%}
\section{Slots}{
  Common to \emph{all} matrix object in the package:
  \describe{
    \item{\code{Dim}:}{Object of class \code{"integer"} - the dimensions
     of the matrix - must be an integer vector with exactly two
     non-negative values.}
    \item{\code{Dimnames}:}{list of length two; each component
      containing NULL or a \code{\link{character}} vector length
      equal the corresponding \code{Dim} element.}
  }
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "matrix", to = "nMatrix")}:
      Note that these coercions (must) coerce \code{\link{NA}}s to
      non-zero, hence conceptually \code{TRUE}.
      This is particularly important when
      \code{\linkS4class{sparseMatrix}} objects are coerced to
      \code{"nMatrix"} and hence to \code{\linkS4class{nsparseMatrix}}.
    }
  }

  --- --- ---

  Additional methods contain group methods, such as
  \describe{
    \item{Ops}{\code{signature(e1 = "nMatrix", e2 = "....")}, \dots}
    \item{Arith}{\code{signature(e1 = "nMatrix", e2 = "....")}, \dots}
    \item{Compare}{\code{signature(e1 = "nMatrix", e2 = "....")}, \dots}
    \item{Logic}{\code{signature(e1 = "nMatrix", e2 = "....")}, \dots}

    \item{Summary}{\code{signature(x = "nMatrix", "....")}, \dots}
  }
}
\seealso{
  The classes \code{\linkS4class{lMatrix}},
  \code{\linkS4class{nsparseMatrix}}, and the mother class,
  \code{\linkS4class{Matrix}}.
}
\examples{
getClass("nMatrix")

L3 <- Matrix(upper.tri(diag(3)))
L3 # an "ltCMatrix"
as(L3, "nMatrix") # -> ntC*

## similar, not using Matrix()
as(upper.tri(diag(3)), "nMatrix")# currently "ngTMatrix"
}
