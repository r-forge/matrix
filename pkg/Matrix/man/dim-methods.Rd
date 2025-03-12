\name{dim-methods}
\title{Dimensions of a Matrix}
%
\docType{methods}
\keyword{array}
\keyword{attribute}
\keyword{methods}
%
\alias{dim}
\alias{dim-methods}
\alias{dim<-}
\alias{dim<--methods}
\alias{length}
\alias{length-methods}
\alias{length<-}
\alias{length<--methods}
%
\alias{dim,Matrix-method}
\alias{dim,MatrixFactorization-method}
\alias{dim,sparseVector-method}
%
\alias{dim<-,denseMatrix,numeric-method}
\alias{dim<-,sparseMatrix,numeric-method}
\alias{dim<-,sparseVector,numeric-method}
%
\alias{length,Matrix-method}
\alias{length,MatrixFactorization-method}
\alias{length,sparseVector-method}
%
\alias{length<-,denseMatrix,numeric-method}
\alias{length<-,sparseMatrix,numeric-method}
\alias{length<-,sparseVector,numeric-method}
%
\description{TODO.}
\usage{
\S4method{dim}{Matrix}(x)
\S4method{dim}{MatrixFactorization}(x)

\S4method{dim}{sparseVector,numeric}(x) <- value

\S4method{length}{Matrix}(x)
\S4method{length}{MatrixFactorization}(x)
\S4method{length}{sparseVector}(x)
}
\arguments{
\item{x}{.}
\item{value}{.}
}
%% \value{}
%% \seealso{}
%% \examples{}

%% \item{dim}{\code{signature(x = "Matrix")}: extract matrix dimensions
%%   \code{\link{dim}}.}

%% \item{dim<-}{\code{signature(x = "Matrix", value = "ANY")}: where
%%   \code{value} is integer of length 2.  Allows to \emph{reshape}
%%   Matrix objects, but only when \code{prod(value) == prod(dim(x))}.}

%% \item{length}{\code{signature(x = "Matrix")}: simply defined as
%%   \code{prod(dim(x))} (and hence of mode \code{"double"}).}

%% \item{dim<-}{\code{signature(x = "sparseMatrix", value = "ANY")}:
%%   allows to \emph{reshape} a sparse matrix to a sparse matrix with
%%   the same entries but different dimensions. \code{value} must be of
%%   length two and fulfill \code{prod(value) == prod(dim(x))}.}

\examples{
(M <- Matrix(c(0, 1, 0, 0), 6, 4))
dim(M)
## can reshape it even :
dim(M) <- c(2, 12)
M
stopifnot(identical(M, Matrix(c(0, 1, 0, 0), 2, 12)))
}
