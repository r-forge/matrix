\name{coerce-methods}
\title{Coercions in Package \pkg{Matrix}}
%
\docType{methods}
\keyword{methods}
%
\alias{coerce}
\alias{coerce-methods}
\alias{as.logical-methods}
\alias{as.integer-methods}
\alias{as.numeric-methods}
\alias{as.complex-methods}
\alias{as.vector-methods}
\alias{as.matrix-methods}
\alias{as.array-methods}
%
\alias{as.logical,Matrix-method}
\alias{as.logical,sparseVector-method}
\alias{as.integer,Matrix-method}
\alias{as.integer,sparseVector-method}
\alias{as.numeric,Matrix-method}
\alias{as.numeric,sparseVector-method}
\alias{as.complex,Matrix-method}
\alias{as.complex,sparseVector-method}
\alias{as.vector,Matrix-method}
\alias{as.vector,sparseVector-method}
\alias{as.matrix,Matrix-method}
\alias{as.matrix,sparseVector-method}
\alias{as.array,Matrix-method}
\alias{as.array,sparseVector-method}
%
\alias{coerce,ANY,Matrix-method}
\alias{coerce,ANY,denseMatrix-method}
\alias{coerce,ANY,sparseMatrix-method}
\alias{coerce,ANY,sparseVector-method}
\alias{coerce,Matrix,CsparseMatrix-method}
\alias{coerce,Matrix,RsparseMatrix-method}
\alias{coerce,Matrix,TsparseMatrix-method}
\alias{coerce,Matrix,copMatrix-method}
\alias{coerce,Matrix,corMatrix-method}
\alias{coerce,Matrix,dMatrix-method}
\alias{coerce,Matrix,ddenseMatrix-method}
\alias{coerce,Matrix,denseMatrix-method}
\alias{coerce,Matrix,diagonalMatrix-method}
\alias{coerce,Matrix,dsparseMatrix-method}
\alias{coerce,Matrix,generalMatrix-method}
\alias{coerce,Matrix,iMatrix-method}
\alias{coerce,Matrix,idenseMatrix-method}
\alias{coerce,Matrix,indMatrix-method}
\alias{coerce,Matrix,isparseMatrix-method}
\alias{coerce,Matrix,lMatrix-method}
\alias{coerce,Matrix,ldenseMatrix-method}
\alias{coerce,Matrix,lsparseMatrix-method}
\alias{coerce,Matrix,matrix-method}
\alias{coerce,Matrix,nMatrix-method}
\alias{coerce,Matrix,ndenseMatrix-method}
\alias{coerce,Matrix,nsparseMatrix-method}
\alias{coerce,Matrix,pMatrix-method}
\alias{coerce,Matrix,packedMatrix-method}
\alias{coerce,Matrix,posdefMatrix-method}
\alias{coerce,Matrix,sparseMatrix-method}
\alias{coerce,Matrix,sparseVector-method}
\alias{coerce,Matrix,symmetricMatrix-method}
\alias{coerce,Matrix,triangularMatrix-method}
\alias{coerce,Matrix,unpackedMatrix-method}
\alias{coerce,Matrix,zMatrix-method}
\alias{coerce,Matrix,zdenseMatrix-method}
\alias{coerce,Matrix,zsparseMatrix-method}
\alias{coerce,diagonalMatrix,symmetricMatrix-method}
\alias{coerce,diagonalMatrix,triangularMatrix-method}
\alias{coerce,dpoMatrix,corMatrix-method}
\alias{coerce,dppMatrix,copMatrix-method}
\alias{coerce,dspMatrix,copMatrix-method}
\alias{coerce,dsyMatrix,corMatrix-method}
\alias{coerce,factor,sparseMatrix-method}
\alias{coerce,generalMatrix,packedMatrix-method}
\alias{coerce,indMatrix,pMatrix-method}
\alias{coerce,indMatrix,packedMatrix-method}
\alias{coerce,list,indMatrix-method}
\alias{coerce,matrix,CsparseMatrix-method}
\alias{coerce,matrix,Matrix-method}
\alias{coerce,matrix,RsparseMatrix-method}
\alias{coerce,matrix,TsparseMatrix-method}
\alias{coerce,matrix,copMatrix-method}
\alias{coerce,matrix,corMatrix-method}
\alias{coerce,matrix,dMatrix-method}
\alias{coerce,matrix,dMatrix-method}
\alias{coerce,matrix,ddenseMatrix-method}
\alias{coerce,matrix,denseMatrix-method}
\alias{coerce,matrix,dgCMatrix-method}
\alias{coerce,matrix,diagonalMatrix-method}
\alias{coerce,matrix,dsparseMatrix-method}
\alias{coerce,matrix,generalMatrix-method}
\alias{coerce,matrix,iMatrix-method}
\alias{coerce,matrix,idenseMatrix-method}
\alias{coerce,matrix,indMatrix-method}
\alias{coerce,matrix,isparseMatrix-method}
\alias{coerce,matrix,lMatrix-method}
\alias{coerce,matrix,lMatrix-method}
\alias{coerce,matrix,ldenseMatrix-method}
\alias{coerce,matrix,lsparseMatrix-method}
\alias{coerce,matrix,nMatrix-method}
\alias{coerce,matrix,ndenseMatrix-method}
\alias{coerce,matrix,nsparseMatrix-method}
\alias{coerce,matrix,pMatrix-method}
\alias{coerce,matrix,packedMatrix-method}
\alias{coerce,matrix,posdefMatrix-method}
\alias{coerce,matrix,sparseMatrix-method}
\alias{coerce,matrix,sparseVector-method}
\alias{coerce,matrix,symmetricMatrix-method}
\alias{coerce,matrix,triangularMatrix-method}
\alias{coerce,matrix,unpackedMatrix-method}
\alias{coerce,matrix,zMatrix-method}
\alias{coerce,matrix,zdenseMatrix-method}
\alias{coerce,matrix,zsparseMatrix-method}
\alias{coerce,nsparseMatrix,indMatrix-method}
\alias{coerce,nsparseMatrix,pMatrix-method}
\alias{coerce,numeric,indMatrix-method}
\alias{coerce,numeric,pMatrix-method}
\alias{coerce,sparseVector,CsparseMatrix-method}
\alias{coerce,sparseVector,Matrix-method}
\alias{coerce,sparseVector,RsparseMatrix-method}
\alias{coerce,sparseVector,TsparseMatrix-method}
\alias{coerce,sparseVector,denseMatrix-method}
\alias{coerce,sparseVector,dsparseVector-method}
\alias{coerce,sparseVector,generalMatrix-method}
\alias{coerce,sparseVector,isparseVector-method}
\alias{coerce,sparseVector,lsparseVector-method}
\alias{coerce,sparseVector,nsparseVector-method}
\alias{coerce,sparseVector,sparseMatrix-method}
\alias{coerce,sparseVector,unpackedMatrix-method}
\alias{coerce,sparseVector,zsparseVector-method}
\alias{coerce,vector,CsparseMatrix-method}
\alias{coerce,vector,Matrix-method}
\alias{coerce,vector,RsparseMatrix-method}
\alias{coerce,vector,TsparseMatrix-method}
\alias{coerce,vector,dMatrix-method}
\alias{coerce,vector,dMatrix-method}
\alias{coerce,vector,ddenseMatrix-method}
\alias{coerce,vector,denseMatrix-method}
\alias{coerce,vector,dsparseMatrix-method}
\alias{coerce,vector,dsparseVector-method}
\alias{coerce,vector,generalMatrix-method}
\alias{coerce,vector,iMatrix-method}
\alias{coerce,vector,idenseMatrix-method}
\alias{coerce,vector,isparseMatrix-method}
\alias{coerce,vector,isparseVector-method}
\alias{coerce,vector,lMatrix-method}
\alias{coerce,vector,lMatrix-method}
\alias{coerce,vector,ldenseMatrix-method}
\alias{coerce,vector,lsparseMatrix-method}
\alias{coerce,vector,lsparseVector-method}
\alias{coerce,vector,nMatrix-method}
\alias{coerce,vector,ndenseMatrix-method}
\alias{coerce,vector,nsparseMatrix-method}
\alias{coerce,vector,nsparseVector-method}
\alias{coerce,vector,sparseMatrix-method}
\alias{coerce,vector,sparseVector-method}
\alias{coerce,vector,unpackedMatrix-method}
\alias{coerce,vector,zMatrix-method}
\alias{coerce,vector,zdenseMatrix-method}
\alias{coerce,vector,zsparseMatrix-method}
\alias{coerce,vector,zsparseVector-method}
%
\description{TODO.}

%% \item{as.matrix, as.array}{\code{signature(x = "Matrix")}: the same as
%%   \code{as(x, "matrix")}; see also the note below.}
%% \item{as.vector}{\code{signature(x = "Matrix", mode = "missing")}:
%%   \code{as.vector(m)} should be identical to \code{as.vector(as(m,
%% 	"matrix"))}, implemented more efficiently for some subclasses.}
%% \item{as(x, "vector"), as(x, "numeric")}{etc, similarly.}
%% \item{coerce}{\code{signature(from = "ANY", to = "Matrix")}: This
%%   relies on a correct \code{\link{as.matrix}()} method for \code{from}.}

%% \item{coerce}{\code{signature(from = "matrix", to = "nMatrix")}:
%%   Note that these coercions (must) coerce \code{\link{NA}}s to
%%   non-zero, hence conceptually \code{TRUE}.
%%   This is particularly important when
%%   \code{\linkS4class{sparseMatrix}} objects are coerced to
%%   \code{"nMatrix"} and hence to \code{\linkS4class{nsparseMatrix}}.}

%% \item{coerce}{\code{signature(from = "factor", to = "sparseMatrix")}:
%%   Coercion of a factor to \code{"sparseMatrix"} produces the matrix
%%   of indicator \bold{rows} stored as an object of class
%%   \code{"dgCMatrix"}.  To obtain columns representing the interaction
%%   of the factor and a numeric covariate, replace the \code{"x"} slot
%%   of the result by the numeric covariate then take the transpose.
%%   Missing values (\code{\link{NA}}) from the factor are translated
%%   to columns of all \code{0}s.}

%% \item{\code{coerce}}{\code{signature(from = "numeric", to = "indMatrix")}:
%%   supporting typical \code{indMatrix} construction from
%%   a vector of positive integers.  Row indexing is assumed.}
%% \item{\code{coerce}}{\code{signature(from = "list", to = "indMatrix")}:
%%   supporting \code{indMatrix} construction for row \emph{and}
%%   column indexing, including index vectors of length 0 and
%%   index vectors whose maximum is less than the number of rows
%%   or columns being indexed.}
%% \item{\code{coerce}}{\code{signature(from = "indMatrix", to = "matrix")}:
%%   coercion to a traditional \code{\link{matrix}} of \link{logical} type,
%%   with \code{FALSE} and \code{TRUE} in place of 0 and 1.}
