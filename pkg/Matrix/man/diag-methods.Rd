\name{diag-methods}
\title{Matrix Diagonals}
%
\docType{methods}
\keyword{algebra}
\keyword{math}
\keyword{methods}
%
\alias{diag}
\alias{diag-methods}
\alias{diag<-}
\alias{diag<--methods}
%
\alias{diag,CsparseMatrix-method}
\alias{diag,RsparseMatrix-method}
\alias{diag,TsparseMatrix-method}
\alias{diag,denseCholesky-method}
\alias{diag,denseMatrix-method}
\alias{diag,diagonalMatrix-method}
\alias{diag,indMatrix-method}
\alias{diag,sparseCholesky-method}
%
\alias{diag<-,CsparseMatrix-method}
\alias{diag<-,RsparseMatrix-method}
\alias{diag<-,TsparseMatrix-method}
\alias{diag<-,denseMatrix-method}
\alias{diag<-,diagonalMatrix-method}
\alias{diag<-,indMatrix-method}
%
\description{TODO.}
\usage{
diag(x = 1, nrow, ncol, names = TRUE)
diag(x) <- value
}
\arguments{
\item{x}{.}
\item{nrow}{.}
\item{ncol}{.}
\item{names}{.}
\item{value}{.}
}
%% \value{}
%% \seealso{}
%% \examples{}

\examples{
(M <- Matrix(c(0, 1, 0, 0), 6, 4))
diag(M)
}
