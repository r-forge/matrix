\name{pMatrix-class}
\docType{class}
\alias{pMatrix-class}
\alias{\%*\%,matrix,pMatrix-method}
\alias{\%*\%,pMatrix,matrix-method}
\alias{coerce,integer,pMatrix-method}
\alias{coerce,pMatrix,matrix-method}
\alias{solve,pMatrix,missing-method}
\alias{t,pMatrix-method}
\title{Permutation matrices}
\description{The \code{"pMatrix"} class is the class of permutation
  matrices stored as 1-based integer permutation vectors.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("pMatrix", ...)}
  or by coercion from an integer permutation vector.
}
\section{Slots}{
  \describe{
    \item{\code{perm}:}{An integer, 1-based permutation vector, i.e.
      an integer vector of length \code{Dim[1]} whose elements form a
      permutation of \code{1:Dim[1]}.}
    \item{\code{Dim}:}{Object of class \code{"integer"}. The dimensions
      of the matrix which must be a two-element vector of equal,
      non-negative integers.}
    \item{\code{Dimnames}:}{list of length two; each component
      containing NULL or a \code{\link{character}} vector length
      equal the corresponding \code{Dim} element.}
  }
}
\section{Extends}{
Class \code{"Matrix"}, directly.
}
\section{Methods}{
  \describe{
    \item{\%*\%}{\code{signature(x = "matrix", y = "pMatrix")}: ... }
    \item{\%*\%}{\code{signature(x = "pMatrix", y = "matrix")}: ... }
    \item{coerce}{\code{signature(from = "integer", to = "pMatrix")}: ... }
    \item{coerce}{\code{signature(from = "pMatrix", to = "matrix")}: ... }
    \item{solve}{\code{signature(a = "pMatrix", b = "missing")}: return
      the inverse permutation matrix.}
    \item{t}{\code{signature(x = "pMatrix")}: return the transpose of
      the permuation matrix (which is also the inverse of the
      permutation matrix).}
  }
}
%\seealso{}
\examples{
(pm1 <- as(as.integer(c(2,3,1)), "pMatrix"))
t(pm1)
solve(pm1)
pm1 \%*\% as(t(pm1), "matrix")   # check that the transpose is the inverse
(mm <- round(array(rnorm(3 * 3), c(3, 3)), 2))
mm \%*\% pm1
pm1 \%*\% mm
try(as(as.integer(c(3,3,1)), "pMatrix"))
}
\keyword{classes}
