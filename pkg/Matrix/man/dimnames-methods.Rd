\name{dimnames-methods}
\title{Dimnames of a Matrix}
%
\docType{methods}
\keyword{array}
\keyword{attribute}
\keyword{methods}
%
\alias{dimnames}
\alias{dimnames-methods}
\alias{dimnames<-}
\alias{dimnames<--methods}
%
\alias{dimnames,Matrix-method}
\alias{dimnames,MatrixFactorization-method}
\alias{dimnames,sparseVector-method}
\alias{dimnames,symmetricMatrix-method}
%
\alias{dimnames<-,Matrix,NULL-method}
\alias{dimnames<-,Matrix,list-method}
\alias{dimnames<-,MatrixFactorization,NULL-method}
\alias{dimnames<-,MatrixFactorization,list-method}
\alias{dimnames<-,generalMatrix,NULL-method}
\alias{dimnames<-,generalMatrix,list-method}
\alias{dimnames<-,symmetricMatrix,NULL-method}
\alias{dimnames<-,symmetricMatrix,list-method}
%
\description{TODO.}
\usage{
\S4method{dimnames}{Matrix}(x)
\S4method{dimnames}{MatrixFactorization}(x)

\S4method{dimnames}{Matrix,list}(x) <- value
\S4method{dimnames}{MatrixFactorization,list}(x) <- value
}
\arguments{
\item{x}{.}
\item{value}{.}
}
%% \value{}
%% \seealso{}
%% \examples{}

%% \item{dimnames}{\code{signature(x = "Matrix")}: extract
%%   \code{\link{dimnames}}.}
%% \item{dimnames<-}{\code{signature(x = "Matrix", value = "list")}: set
%%   the \code{dimnames} to a \code{\link{list}} of length 2, see
%%   \code{\link{dimnames<-}}.}
%% \item{dimnames}{\code{signature(object = "symmetricMatrix")}:
%%    returns \emph{symmetric} \code{\link{dimnames}}, even when the
%%    \code{Dimnames} slot only has row or column names.  This allows to
%%    save storage for large (typically sparse) symmetric matrices.}
