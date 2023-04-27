\name{Schur-class}
\title{Schur Factorizations}
%
\docType{class}
\keyword{algebra}
\keyword{array}
\keyword{classes}
%
\alias{Schur-class}
%
\alias{show,Schur-method}
%
\description{
  \code{Schur} is the class of Schur factorizations of real, square
  matrices.  The Schur factorization generalizes the eigendecomposition
  of diagonalizable (i.e., non-defective) square matrices to general
  square matrices.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}}{inherited from virtual class
      \code{\linkS4class{MatrixFactorization}}.}
    \item{\code{Q}}{an orthogonal matrix,
      inheriting from virtual class \code{\linkS4class{Matrix}}.}
    \item{\code{T}}{a block upper triangular matrix,
      inheriting from virtual class \code{\linkS4class{Matrix}}.
      The diagonal blocks have dimensions 1-by-1 or 2-by-2.}
    \item{\code{EValues}}{a numeric or complex vector containing
      the eigenvalues of the diagonal blocks of \code{T}, which are
      the eigenvalues of \code{T} and consequently of the factorized
      matrix.}
  }
}
\section{Extends}{
  Class \code{\linkS4class{MatrixFactorization}}, directly.
}
\section{Instantiation}{
  Objects can be generated directly by calls of the form
  \code{new("Schur", ...)}, but they are more typically obtained
  as the value of \code{\link{Schur}(x)} for \code{x} inheriting from
  \code{\linkS4class{Matrix}} (often \code{\linkS4class{dgeMatrix}}).
}
\section{Methods}{
  \describe{
    \item{expand2}{\code{signature(x = "Schur")}:
      see \code{\link{expand2-methods}}.}
    \item{show}{\code{signature(object = "Schur")}:
      prints the internal representation of the factorization using
      \code{\link{str}}.}
  }
}
\seealso{
  Class \code{\linkS4class{dgeMatrix}}.
  
  Generic functions \code{\link{expand2}} and \code{\link{Schur}}.
}
\examples{
showClass("Schur")
Schur(M <- Matrix(c(1:7, 10:2), 4,4))
## Trivial, of course:
str(Schur(Diagonal(5)))

## for more examples, see help("Schur")
}
