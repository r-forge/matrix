\name{lgCMatrix-class}
\docType{class}
\alias{lgCMatrix-class}
\alias{ltCMatrix-class}
\alias{lsCMatrix-class}
\alias{\%*\%,lgCMatrix,lgCMatrix-method}
\alias{coerce,matrix,lgCMatrix-method}
\alias{t,lgCMatrix-method}
\title{Compressed, sparse, column-oriented logical matrices}
\description{The \code{lgCMatrix} class is a class of general, sparse,
  logical matrices in the compressed, sparse, column-oriented format.
  \code{lsCMatrix} and \code{ltCMatrix} are similar classes representing
  symmetric and triangular matrices respectively.

  In this implementation the non-zero elements in the columns are sorted into
  increasing row order.  There is no \code{x} slot in this class because
  the only non-zero value stored is \code{TRUE}.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("lgCMatrix", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{uplo}:}{Object of class \code{"character"}. Must be
      either "U", for upper triangular, and "L", for lower
      triangular. Present in the triangular and symmetric classes but not
      in the general class.}
    \item{\code{diag}:}{Object of class \code{"character"}. Must be
      either \code{"U"}, for unit triangular (diagonal is all ones), or
      \code{"N"} for non-unit.  The implicit diagonal elements are not
      explicitly stored when \code{diag} is \code{"U"}.  Present in the
      triangular classes only.}
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
  \code{\link{dgCMatrix-class}}
}
%\examples{}
\keyword{classes}
\keyword{algebra}
