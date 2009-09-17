\name{chol2inv-methods}
\docType{methods}
\alias{chol2inv-methods}
\alias{chol2inv,ANY-method}
\alias{chol2inv,dtrMatrix-method}
\title{Inverse from Choleski or QR Decomposition -- Matrix Methods}
\description{
  Invert a symmetric, positive definite square matrix from its Choleski
  decomposition.  Equivalently, compute \eqn{(X'X)^{-1}}{(X'X)^(-1)}
  from the (\eqn{R} part) of the QR decomposition of \eqn{X}.
  \cr
  Even more generally, given an upper triangular matrix \eqn{R},
  compute \eqn{(R'R)^{-1}}{(R'R)^(-1)}.
}
\usage{
chol2inv(x)
}
\section{Methods}{
  \describe{
    \item{x = "ANY"}{the default method from \pkg{base}, see
      \code{\link[base]{chol2inv}}, for traditional matrices.}

    \item{x = "dtrMatrix"}{method for the numeric triangular matrices,
      built on the same LAPACK \command{DPOTRI} function as the base
      method.}
  }
}
\seealso{
  \code{\link{chol}} (for \code{\linkS4class{Matrix}} objects),
  \code{\link{solve}}.
}
\examples{
(M  <- Matrix(cbind(1, 1:3, c(1,3,7))))
(cM <- chol(M)) # a "Cholesky" object, inheriting from "dtrMatrix"
chol2inv(cM) \%*\% M # the identity
stopifnot(all(chol2inv(cM) \%*\% M - Diagonal(nrow(M))) < 1e-10)
}
\keyword{methods}
\keyword{algebra}

