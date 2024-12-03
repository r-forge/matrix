\name{RsparseMatrix-class}
\title{Class "RsparseMatrix" of Sparse Matrices in Row-compressed Form}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{RsparseMatrix-class}
%
\alias{ngRMatrix-class}
\alias{lgRMatrix-class}
\alias{igRMatrix-class}
\alias{dgRMatrix-class}
\alias{zgRMatrix-class}
%
\alias{nsRMatrix-class}
\alias{lsRMatrix-class}
\alias{isRMatrix-class}
\alias{dsRMatrix-class}
\alias{zsRMatrix-class}
%
\alias{ntRMatrix-class}
\alias{ltRMatrix-class}
\alias{itRMatrix-class}
\alias{dtRMatrix-class}
\alias{ztRMatrix-class}
%
\description{The \code{"RsparseMatrix"} class is the virtual class of
  all sparse matrices coded in sorted compressed row-oriented form.
  Since it is a virtual class, no objects may be created from it.  See
  \code{showClass("RsparseMatrix")} for its subclasses.
}
\section{Slots}{
  \describe{
    \item{\code{j}:}{Object of class \code{"integer"} of length \code{nnzero}
      (number of non-zero elements).  These are the row numbers for
      each non-zero element in the matrix.}
    \item{\code{p}:}{Object of class \code{"integer"} of pointers, one
      for each row, to the initial (zero-based) index of elements in
      the row.}
    \item{\code{Dim}, \code{Dimnames}:}{inherited from
      the superclass, see \code{\linkS4class{sparseMatrix}}.}
  }
}
\section{Extends}{
  Class \code{"sparseMatrix"}, directly.
  Class \code{"Matrix"}, by class \code{"sparseMatrix"}.
}
\section{Methods}{
  Originally, \bold{few} methods were defined on purpose, as we
  rather use the \code{\linkS4class{CsparseMatrix}} in \pkg{Matrix}.
  Then, more methods were added but \emph{beware} that these
  typically do \emph{not} return \code{"RsparseMatrix"} results, but
  rather Csparse* or Tsparse* ones; e.g., \code{R[i, j] <- v} for an
  \code{"RsparseMatrix"} \code{R} works, but after the assignment, \code{R}
  is a (triplet) \code{"TsparseMatrix"}.
  \describe{
    \item{t}{\code{signature(x = "RsparseMatrix")}: ... }
    \item{coerce}{\code{signature(from = "RsparseMatrix", to = "CsparseMatrix")}: ... }
    \item{coerce}{\code{signature(from = "RsparseMatrix", to = "TsparseMatrix")}: ... }
  }
}
\seealso{
  its superclass, \code{\linkS4class{sparseMatrix}}, and, e.g., class
  \code{\linkS4class{dgRMatrix}} for the links to other classes.
}
\examples{
showClass("RsparseMatrix")

%% moved from ./dsRMatrix-class.Rd :
(m0 <- new("dsRMatrix"))
m2 <- new("dsRMatrix", Dim = c(2L,2L),
          x = c(3,1), j = c(1L,1L), p = 0:2)
m2
stopifnot(colSums(as(m2, "TsparseMatrix")) == 3:4)
str(m2)
(ds2 <- forceSymmetric(diag(2))) # dsy*
dR <- as(ds2, "RsparseMatrix")
dR # dsRMatrix

%% moved from ./dtRMatrix-class.Rd :
(m0 <- new("dtRMatrix"))
(m2 <- new("dtRMatrix", Dim = c(2L,2L),
                        x = c(5, 1:2), p = c(0L,2:3), j= c(0:1,1L)))
str(m2)
(m3 <- as(Diagonal(2), "RsparseMatrix"))# --> dtRMatrix
}
