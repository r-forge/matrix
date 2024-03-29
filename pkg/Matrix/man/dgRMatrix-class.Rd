\name{dgRMatrix-class}
\title{Sparse Compressed, Row-oriented Numeric Matrices}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{dgRMatrix-class}
%
\alias{determinant,dgRMatrix,logical-method}
%
\description{The \code{dgRMatrix} class is a class of sparse numeric
  matrices in the compressed, sparse, row-oriented format.  In this
  implementation the non-zero elements in the rows are sorted into
  increasing column order.

  \bold{Note:} The column-oriented sparse classes, e.g.,
  \code{\linkS4class{dgCMatrix}}, are preferred and better supported in
  the \pkg{Matrix} package.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("dgRMatrix", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{j}:}{Object of class \code{"integer"} of length nnzero
      (number of non-zero elements).  These are the column numbers for
      each non-zero element in the matrix.}
    \item{\code{p}:}{Object of class \code{"integer"} of pointers, one
      for each row, to the initial (zero-based) index of elements in
      the row.}
    \item{\code{x}:}{Object of class \code{"numeric"} - the non-zero
      elements of the matrix.}
    \item{\code{Dim}:}{Object of class \code{"integer"} - the dimensions
      of the matrix.}
  }
}
\section{Methods}{
  \describe{
    \item{diag}{\code{signature(x = "dgRMatrix")}: returns the diagonal
      of \code{x}}
    \item{dim}{\code{signature(x = "dgRMatrix")}: returns the dimensions
      of \code{x}}
    \item{image}{\code{signature(x = "dgRMatrix")}: plots an image of
      \code{x} using the \code{\link[lattice]{levelplot}} function}
  }
}
\seealso{
  the \code{\linkS4class{RsparseMatrix}} class, the virtual class of all
  sparse compressed \bold{r}ow-oriented matrices, with its methods.
  The \code{\linkS4class{dgCMatrix}} class (\bold{c}olumn compressed
  sparse) is really preferred.
}
