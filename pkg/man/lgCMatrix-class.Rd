\name{lgCMatrix-class}
\docType{class}
\alias{lgCMatrix-class}
\alias{\%*\%,lgCMatrix,lgCMatrix-method}
\alias{coerce,matrix,lgCMatrix-method}
\alias{t,lgCMatrix-method}
\title{Compressed, sparse, column-oriented logical matrices}
\description{The \code{dgCMatrix} class is a class of sparse logical
  matrices in the compressed, sparse, column-oriented format.  In this
  implementation the non-zero elements in the columns are sorted into
  increasing row order.  There is no \code{x} slot in this class because
  the only non-zero value stored is \code{TRUE}.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("lgCMatrix", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{p}:}{Object of class \code{"integer"} of pointers, one
      for each column, to the initial (zero-based) index of elements in
      the column.}
    \item{\code{i}:}{Object of class \code{"integer"} of length nnzero
      (number of non-zero elements).  These are the row numbers for
      each TRUE element in the matrix.  All other elements are FALSE.}
    \item{\code{Dim}:}{Object of class \code{"integer"} - the dimensions
      of the matrix.}
  }
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "dgCMatrix", to = "lgCMatrix")}}
    \item{t}{\code{signature(x = "lgCMatrix")}: returns the transpose
      of \code{x}}
  }
}
%\references{}
%\author{}
\note{
  This class is most frequently used during the symbolic analysis phase
  of an algorithm involving sparse matrices.  Frequently such algorithms
  involve two phases: a symbolic phase wherein the positions of the
  non-zeros in the result are determined and a numeric phase wherein the
  actual results are calculated.  During the symbolic phase only the
  positions of the non-zero elements in any operands are of interest,
  hence any numeric sparse matrices can be treated as a logical sparse
  matrices.
}
\seealso{
  \code{\link{dgCMatrix-class}}, \code{\link{ltCMatrix-class}},
    \code{\link{lsCMatrix-class}}
}
%\examples{}
\keyword{classes}
\keyword{algebra}
