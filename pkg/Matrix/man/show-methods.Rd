\name{show-methods}
\title{Show a Matrix}
%
\docType{methods}
\keyword{character}
\keyword{print}
\keyword{methods}
%
\alias{show}
\alias{show-methods}
%
\alias{show,BunchKaufman-method}
\alias{show,Cholesky-method}
\alias{show,LU-method}
\alias{show,MatrixFactorization-method}
\alias{show,QR-method}
\alias{show,Schur-method}
\alias{show,denseMatrix-method}
\alias{show,diagonalMatrix-method}
\alias{show,sparseMatrix-method}
\alias{show,sparseVector-method}
%
\alias{print,diagonalMatrix-method}
\alias{print,sparseMatrix-method}
%
\alias{format,sparseMatrix-method}
%
\alias{summary,diagonalMatrix-method}
\alias{summary,sparseMatrix-method}
%
\alias{print.diagSummary}
\alias{print.sparseSummary}
%
\description{TODO.}
\usage{
show(object)
}
\arguments{
\item{object}{.}
}
%% \value{}
%% \seealso{}
%% \examples{}

%% \item{show}{\code{signature(object = "Matrix")}: \code{\link{show}}
%%   method for \code{\link{print}}ing.  For printing \emph{sparse}
%%   matrices, see \code{\link{printSpMatrix}}.}

%% \item{show}{\code{(object = "sparseMatrix")}: The
%%   \code{\link{show}} method for sparse matrices prints
%%   \emph{\dQuote{structural}} zeroes as \code{"."} using
%%   \code{\link{printSpMatrix}()} which allows further customization.}

%% \item{print}{\code{signature(x = "sparseMatrix")}, ....\cr
%%   The \code{\link{print}} method for sparse matrices by default is the
%%   same as \code{show()} but can be called with extra optional
%%   arguments, see \code{\link{printSpMatrix}()}.}

%% \item{format}{\code{signature(x = "sparseMatrix")}, ....\cr
%%   The \code{\link{format}} method for sparse matrices, see
%%   \code{\link{formatSpMatrix}()} for details such as the extra
%%   optional arguments.}

%% \item{summary}{\code{(object = "sparseMatrix", uniqT=FALSE)}: Returns
%%   an object of S3 class \code{"sparseSummary"} which is basically a
%%   \code{\link{data.frame}} with columns \code{(i,j,x)} (or just
%%   \code{(i,j)} for \code{\linkS4class{nsparseMatrix}} class objects)
%%   with the stored (typically non-zero) entries.  The
%%   \code{\link{print}} method resembles Matlab's way of printing
%%   sparse matrices, and also the MatrixMarket format, see
%%   \code{\link{writeMM}}.}
