\name{dtCMatrix-class}
\docType{class}
\alias{dtCMatrix-class}
\alias{t,dtCMatrix-method}
\alias{coerce,dtCMatrix,dgTMatrix-method}
\alias{coerce,dtCMatrix,dgeMatrix-method}
\title{Triangular, compressed, sparse column matrices}
\description{The \code{dsCMatrix} class is a class of triangular, sparse
  matrices in the compressed, column-oriented format.  In this
  implementation the non-zero elements in the columns are sorted into
  increasing row order.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("dtCMatrix", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{uplo}:}{A character object indicating if the matrix is upper
      triangular (\code{"U"} or \code{"u"}) or lower triangular
      (\code{"L"} or \code{"l"}).}
    \item{\code{p}:}{Object of class \code{"integer"} of pointers, one
      for each column, to the initial (zero-based) index of elements in
      the column.}
    \item{\code{i}:}{Object of class \code{"integer"} of length nnzero
      (number of non-zero elements).  These are the row numbers for
      each non-zero element in the matrix.}
    \item{\code{x}:}{Object of class \code{"numeric"} - the non-zero
      elements of the matrix.}
    \item{\code{factors}:}{Object of class \code{"list"} - a list
      of factorizations of the matrix. }
    \item{\code{Dim}:}{Object of class \code{"integer"} - the dimensions
     of the matrix - must be an integer vector with exactly two
     non-negative values.}
  }
}
\section{Extends}{
  Class \code{"dgCMatrix"}, directly.
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "dtCMatrix", to = "dgTMatrix")}}
    \item{coerce}{\code{signature(from = "dtCMatrix", to = "dgeMatrix")}}
    \item{t}{\code{signature(x = "dtCMatrix")}: returns the transpose of
      \code{x}}
  }
}
%\references{}
%\author{}
%\note{}
\seealso{
  \code{\link{dgCMatrix-class}}, \code{\link{dgTMatrix-class}},
  \code{\link{dgeMatrix-class}}
}
%\examples{}
\keyword{classes}
\keyword{algebra}
