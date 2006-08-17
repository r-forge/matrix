\name{sparseQR-class}
\docType{class}
\alias{sparseQR-class}
\title{Sparse QR decomposition of a sparse matrix}
\description{Objects of this class represent a QR decomposition of a
  sparse rectangular matrix.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("sparseQR", ...)}
  but are more commonly created by function \code{\link[base]{qr}} applied
  to a sparse matrix such as a matrix of class
  \code{\linkS4class{dgCMatrix}}.}
}
\section{Slots}{
  \describe{
    \item{\code{V}:}{Object of class \code{"dgCMatrix"} The columns of
      \code{V} are the vectors that generate the Householder
      transformations of which the matrix Q is composed.}
    \item{\code{beta}:}{Object of class \code{"numeric"} The normalizing
      factors for the Householder transformations.}
    \item{\code{p}:}{Object of class \code{"integer"} Permutation
      applied to the rows of the original matrix.}
    \item{\code{R}:}{Object of class \code{"dgCMatrix"} An upper
      triangular matrix of dimension \ }
    \item{\code{q}:}{Object of class \code{"integer"} Permutation
      applied from the right.}
  }
}
\section{Methods}{
No methods defined with class "sparseQR" in the signature.
}
%\references{}
%\author{}
\note{The decomposition is of the form
  \deqn{A = PQRS}{A = PQRS}
  where A is a sparse \eqn{m\times n}{m by n} matrix (\eqn{m\ge n}{m >= n}),
  \eqn{P} is an \eqn{m\times m}{m by m} permutation matrix, \eqn{Q} is a
  \eqn{m\times m}{m by m} orthogonal matrix, \eqn{R} is an \eqn{m\times
    n}{m by n} matrix that is zero below the main diagonal and \eqn{S}
  is a \eqn{n\times n}{n by n} permutation matrix.  The permutation
  matrices are chosen to reduce the amount of "fill-in" in the matrix
  \eqn{R}.

  The matrix \eqn{Q} is a "virtual matrix".  It is the product of
  \eqn{n} Householder transformations.  The information to generate
  these Householder transformations is stored in the \code{V} and
  \code{beta} slots.
}
\seealso{
  \code{\link[base]{qr}}, \code{\link[base]{qr.Q}},
  \code{\link[base]{qr.R}}, \code{\link[base]{qr.fitted}},
  \code{\link[base]{qr.resid}}, \code{\link[base]{qr.coef}},
  \code{\link[base]{qr.qty}}, \code{\link[base]{qr.qy}},
  \code{\linkS4class{dgCMatrix}}.
}
\examples{
data(KNex); mm <- KNex$mm
str(mmQR <- qr(mm))
}
\keyword{classes}
\keyword{algebra}
\keyword{array}
