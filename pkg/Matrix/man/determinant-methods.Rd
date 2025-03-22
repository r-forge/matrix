\name{determinant-methods}
\title{Matrix Determinant}
%
\docType{methods}
\keyword{algebra}
\keyword{math}
\keyword{methods}
%
\alias{determinant}
\alias{determinant-methods}
%
\alias{determinant,Matrix,missing-method}
\alias{determinant,MatrixFactorization,missing-method}
\alias{determinant,denseBunchKaufman,logical-method}
\alias{determinant,denseCholesky,logical-method}
\alias{determinant,denseLU,logical-method}
\alias{determinant,denseSchur,logical-method}
\alias{determinant,diagonalMatrix,logical-method}
\alias{determinant,generalMatrix,logical-method}
\alias{determinant,indMatrix,logical-method}
\alias{determinant,pMatrix,logical-method}
\alias{determinant,sparseCholesky,logical-method}
\alias{determinant,sparseLU,logical-method}
\alias{determinant,sparseQR,logical-method}
\alias{determinant,symmetricMatrix,logical-method}
\alias{determinant,triangularMatrix,logical-method}
%
\alias{det}
%
\description{TODO.}
\usage{
det(x, \dots)
determinant(x, logarithm = TRUE, \dots)
}
\arguments{
\item{x}{.}
\item{logarithm}{.}
\item{\dots}{.}
}
%% \value{}
%% \seealso{}
%% \examples{}

%% \item{determinant}{\code{signature(x = "Matrix", logarithm = "missing")}: and}
%% \item{determinant}{\code{signature(x = "Matrix", logarithm = "logical")}:
%%   compute the (\eqn{\log}) determinant of \code{x}.  The method
%%   chosen depends on the actual Matrix class of \code{x}.  Note that
%%   \code{\link[base]{det}} also works for all our matrices, calling the
%%   appropriate \code{determinant()} method.  The \code{Matrix::det}
%%   is an exact copy of \code{base::det}, but in the correct
%%   namespace, and hence calling the S4-aware version of
%%   \code{determinant()}.).}

%% \item{determinant}{\code{(x = "sparseMatrix", logarithm=TRUE)}:
%%   \code{\link{determinant}()} methods for sparse matrices typically
%%   work via \code{\link{Cholesky}} or \code{\link{lu}} decompositions.}
%% \item{determinant}{\code{signature(x = "dpoMatrix")}:
%%   Returns the \code{\link{determinant}} of \code{x}, via
%%   \code{chol(x)}, see above.}

\examples{
(M <- Matrix(c(0, 1, 0, 0), 6, 4))
cm <- M[1:4, ] + 10 * Diagonal(4)
stopifnot(all.equal(det(cm),
                    determinant(as(cm, "matrix"), log = FALSE)$modulus,
                    check.attributes = FALSE))
}
