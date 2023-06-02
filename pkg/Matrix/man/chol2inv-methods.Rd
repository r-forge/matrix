\name{chol2inv-methods}
\title{Inverse from Cholesky (or QR) Decomposition}
%
\docType{methods}
\keyword{algebra}
\keyword{array}
\keyword{methods}
%
\alias{chol2inv}
\alias{chol2inv-methods}
%
\alias{chol2inv,ANY-method}
\alias{chol2inv,ddiMatrix-method}
\alias{chol2inv,denseMatrix-method}
\alias{chol2inv,dtCMatrix-method}
\alias{chol2inv,dtpMatrix-method}
\alias{chol2inv,dtrMatrix-method}
\alias{chol2inv,sparseMatrix-method}
%
\description{
  For \emph{formally} upper and lower triangular matrices
  \eqn{U} and \eqn{L}, compute \eqn{(U' U)^{-1}}{(U' U)^(-1)}
  and \eqn{(L L')^{-1}}{(L L')^(-1)}, respectively.
  For all other matrices, consider only the upper triangular
  part \eqn{U} and compute \eqn{(U' U)^{-1}}{(U' U)^(-1)}.
  
  This function can be seen as way to invert a symmetric,
  positive definite matrix from its Cholesky factor.
  Equivalently, it can be seen as a way to compute
  \eqn{(X' X)^{-1}}{(X' X)^(-1)} from the \eqn{R} part
  of the QR factorization of \eqn{X}.
}
% \usage{
% chol2inv(x, \dots)
% }
% \arguments{
%   \item{x}{a matrix(-like) object; see below.}
%   \item{\dots}{not used here; for compatibility with other methods.}
% }
\section{Methods}{
  \describe{
    \item{x = "ANY"}{the default method from \pkg{base}, see
      \code{\link[base]{chol2inv}}, for traditional matrices.}

    \item{x = "dtrMatrix"}{method for the numeric triangular matrices,
      built on the same LAPACK \command{DPOTRI} function as the base
      method.}

    \item{x = "denseMatrix"}{if \code{x} is coercable to a
      \code{\linkS4class{triangularMatrix}}, call the \code{"dtrMatrix"}
      method above.}

    \item{x = "sparseMatrix"}{if \code{x} is coercable to a
      \code{\linkS4class{triangularMatrix}}, use \code{\link{solve}()}
      currently.} %% better algorithms are welcome!
  }
}
\seealso{
  \code{\link{chol}} (for \code{\linkS4class{Matrix}} objects);
  further, \code{\link[base]{chol2inv}} (from the \pkg{base} package),
  \code{\link{solve}}.
}
\examples{
(M  <- Matrix(cbind(1, 1:3, c(1,3,7))))
(cM <- chol(M)) # a "Cholesky" object, inheriting from "dtrMatrix"
chol2inv(cM) \%*\% M # the identity
stopifnot(all(chol2inv(cM) \%*\% M - Diagonal(nrow(M))) < 1e-10)
}
