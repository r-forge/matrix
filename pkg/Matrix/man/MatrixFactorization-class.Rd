\name{MatrixFactorization-class}
\title{Virtual Class "MatrixFactorization" of Matrix Factorizations}
%
\docType{class}
\keyword{algebra}
\keyword{array}
\keyword{classes}
%
\alias{MatrixFactorization-class}
%
\alias{BunchKaufmanFactorization-class}
\alias{CholeskyFactorization-class}
\alias{SchurFactorization-class}
\alias{LU-class}
\alias{QR-class}
%
\alias{determinant,MatrixFactorization,missing-method}
\alias{dim,MatrixFactorization-method}
\alias{dimnames,MatrixFactorization-method}
\alias{show,MatrixFactorization-method}
%
\alias{show,BunchKaufmanFactorization-method}
\alias{show,CholeskyFactorization-method}
\alias{show,SchurFactorization-method}
\alias{show,LU-method}
\alias{show,QR-method}
%
\description{
  \code{MatrixFactorization} is the virtual class of matrix factorizations.
  
  \code{CholeskyFactorization}, \code{BunchKaufmanFactorization},
  \code{SchurFactorization}, \code{LU}, and \code{QR} are the virtual
  subclasses of \code{MatrixFactorization} containing all Cholesky,
  Bunch-Kaufman, Schur, LU, and QR factorizations, respectively.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}}{an integer vector of length 2 giving the
      dimensions of the factorized matrix.}
    \item{\code{Dimnames}}{a list of length 2 preserving the
      \code{dimnames} of the factorized matrix.  Each element
      must be \code{NULL} or a character vector of length equal
      to the corresponding element of \code{Dim}.}
  }
}
\section{Methods}{
  \describe{
    \item{determinant}{\code{signature(x = "MatrixFactorization", logarithm = "missing")}:
      sets \code{logarithm = TRUE} and recalls the generic function.}
    \item{dim}{\code{signature(x = "MatrixFactorization")}:
      returns \code{x@Dim}.}
    \item{dimnames}{\code{signature(x = "MatrixFactorization")}:
      returns \code{x@Dimnames}.}
    \item{show}{\code{signature(object = "MatrixFactorization")}:
      prints the internal representation of the factorization using
      \code{\link{str}}.}
    \item{solve}{\code{signature(a = "MatrixFactorization", b = .)}:
      see \code{\link{solve-methods}}.}
  }
}
\seealso{
  The virtual class \code{\linkS4class{compMatrix}}
  of \emph{factorizable} matrices.
  
  Classes extending \code{CholeskyFactorization}, namely
  \code{\linkS4class{Cholesky}}, \code{\linkS4class{pCholesky}},
  and \code{\linkS4class{CHMfactor}}.

  Classes extending \code{BunchKaufmanFactorization}, namely
  \code{\linkS4class{BunchKaufman}} and \code{\linkS4class{pBunchKaufman}}.

  Classes extending \code{SchurFactorization}, namely
  \code{\linkS4class{Schur}}.

  Classes extending \code{LU}, namely
  \code{\linkS4class{denseLU}} and \code{\linkS4class{sparseLU}}.
  
  Classes extending \code{QR}, namely \code{\linkS4class{sparseQR}}.
  
  Generic functions \code{\link{chol}}, \code{\link{Cholesky}},
  \code{\link{BunchKaufman}}, \code{\link{Schur}}, \code{\link{lu}},
  and \code{\link{qr}} for \emph{computing} factorizations.

  Generic function \code{\link{expand2}} for constructing
  matrix factors from \code{MatrixFactorization} objects.
}
\examples{
showClass("MatrixFactorization")
}
