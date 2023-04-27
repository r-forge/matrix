\name{CHMfactor-class}
\title{Sparse Cholesky Factorizations}
%
\docType{class}
\keyword{algebra}
\keyword{array}
\keyword{classes}
\keyword{programming}
\keyword{utilities}
%
\alias{CHMfactor-class}
\alias{CHMsimpl-class}
\alias{CHMsuper-class}
\alias{dCHMsimpl-class}
\alias{dCHMsuper-class}
\alias{nCHMsimpl-class}
\alias{nCHMsuper-class}
%
\alias{coerce,CHMfactor,CsparseMatrix-method}
\alias{coerce,CHMfactor,Matrix-method}
\alias{coerce,CHMfactor,RsparseMatrix-method}
\alias{coerce,CHMfactor,TsparseMatrix-method}
\alias{coerce,CHMfactor,dMatrix-method}
\alias{coerce,CHMfactor,dsparseMatrix-method}
\alias{coerce,CHMfactor,pMatrix-method}
\alias{coerce,CHMfactor,sparseMatrix-method}
\alias{coerce,CHMfactor,triangularMatrix-method}
\alias{determinant,CHMfactor,logical-method}
\alias{update,CHMfactor-method}
%
\alias{isLDL}
%
\description{
  \code{CHMfactor} is the virtual class of sparse Cholesky
  factorizations of real, symmetric, positive semidefinite matrices.
  The implementation is based on CHOLMOD's \code{cholmod_factor_struct}.
  
  Virtual subclasses \code{CHMsimpl} and \code{CHMsuper} separate the
  simplicial and supernodal variants.  These have nonvirtual subclasses
  \code{[dn]CHMsimpl} and \code{[dn]CHMsuper},
  where prefix \samp{d} and prefix \samp{n} are reserved for numeric and
  symbolic factorizations, respectively.
}
\usage{
isLDL(x)
}
\arguments{
  \item{x}{a Cholesky factorization inheriting from virtual class
    \code{CHMfactor}, almost always the result of a call to generic
    function \code{\link{Cholesky}}.}
}
\section{Slots}{
  Of \code{CHMfactor}:
  \describe{
    \item{\code{Dim}, \code{Dimnames}}{inherited from virtual class
      \code{\linkS4class{MatrixFactorization}}.}
    \item{\code{colcount}}{an integer vector of length \code{Dim[1]}
      giving an \emph{estimate} of the number of nonzero entries in
      each column of the lower triangular Cholesky factor.
      If symbolic analysis was performed prior to factorization,
      then the estimate is exact.}
    \item{\code{perm}}{a integer vector specifying the permutation
      applied to the rows and columns of the factorized matrix.
      If symbolic analysis was performed prior to factorization,
      then the permutation is chosen to reduce fill-in.}
    \item{\code{type}}{an integer vector specifying the factorization
      type.  The details of this vector are subject to change.
      Currently, the length is at least 4 and the first 4 elements
      correspond to members \code{ordering}, \code{is_ll},
      \code{is_super}, and \code{is_monotonic} of the original
      \code{cholmod_factor_struct}.  Subclasses may employ additional
      elements.}
  }
  Of \code{CHMsimpl} (all unused by \code{nCHMsimpl}):
  \describe{
    \item{\code{nz}}{an integer vector of length \code{Dim[1]}
      giving the number of nonzero entries in each column of the
      lower triangular Cholesky factor.  There is at least one
      nonzero entry in each column, because the diagonal elements
      of the factor are stored explicitly.}
    \item{\code{p}}{an integer vector of length \code{Dim[1]+1}.
      Row indices of nonzero entries in column \code{j} of the
      lower triangular Cholesky factor are obtained as
      \code{i[p[j]+seq_len(nz[j])]+1}.}
    %% MJ: Does CHOLMOD document the expected value of p[length(p)] ??
    \item{\code{i}}{an integer vector of length greater than or equal
      to \code{sum(nz)} containing the row indices of nonzero entries
      in the lower triangular Cholesky factor.  These are grouped by
      column and sorted within columns, but the columns themselves
      need not be ordered monotonically.  Columns may be overallocated,
      i.e., the number of elements of \code{i} reserved for column
      \code{j} may exceed \code{nz[j]}.}
  }
  Of \code{dCHMsimpl}:
  \describe{
    \item{\code{x}}{a numeric vector parallel to \code{i} (and of
      equal length) containing the corresponding nonzero entries
      of the lower triangular Cholesky factor.}
  }
  Of \code{CHMsuper}:
  \describe{
    \item{\code{super}, \code{pi}, \code{px}}{integer vectors of
      length \code{nsuper+1}, where \code{nsuper} is the number of
      supernodes.  \code{super[j]+1} is the index of the leftmost
      column of supernode \code{j}.  The row indices of supernode
      \code{j} are obtained as \code{s[pi[j]+seq_len(pi[j+1]-pi[j])]+1}.
      The numeric entries of supernode \code{j} are obtained as
      \code{x[px[j]+seq_len(px[j+1]-px[j])]+1} (if slot \code{x}
      is available).}
    \item{\code{s}}{an integer vector of length greater than or equal
      to \code{Dim[1]} containing the row indices of the supernodes.
      \code{s} may contain duplicates, but not within a supernode,
      where the row indices are strictly increasing.}
  }
  Of \code{dCHMsuper}:
  \describe{
    \item{\code{x}}{a numeric vector of length less than or equal to
      \code{prod(Dim)} containing the numeric entries of the supernodes.}
  }
}
\section{Extends}{
  Class \code{\linkS4class{MatrixFactorization}}, directly.
}
\section{Instantiation}{
  Objects can be generated directly by calls of the form
  \code{new("dCHMsimpl", ...)}, etc., but \code{dCHMsimpl} and
  \code{dCHMsuper} are more typically obtained as the value of
  \code{\link{Cholesky}(x, ...)} for \code{x} inheriting from
  \code{\linkS4class{sparseMatrix}}
  (often \code{\linkS4class{dsCMatrix}}).
  
  There is currently no API outside of calls to \code{\link{new}}
  for generating \code{nCHMsimpl} and \code{nCHMsuper}.
}
\section{Methods}{
  \describe{
    \item{chol2inv}{\code{signature(x = "CHMfactor")}:
      see \code{\link{chol2inv-methods}}.}
    \item{coerce}{\code{signature(from = "CHMfactor", to = "Matrix")}:
      returns a \code{\linkS4class{dtCMatrix}} representing
      the nonunit lower triangular Cholesky factor \eqn{L} in
      \eqn{A = P' L L' P}, \emph{even if} \code{isLDL(from)}
      is \code{TRUE}.}
    \item{coerce}{\code{signature(from = "CHMfactor", to = "pMatrix")}:
      returns a \code{\linkS4class{pMatrix}} representing
      the permutation matrix \eqn{P} in \eqn{A = P' L L' P}.}
    \item{determinant}{\code{signature(x = "CHMfactor", logarithm = "logical")}:
      returns the determinant of the nonunit lower triangular
      Cholesky factor \eqn{L} in \eqn{A = P' L L' P},
      which is the square root of the determinant of \eqn{A},
      \emph{even if} \code{isLDL(x)} is \code{TRUE}.}
    \item{expand}{\code{signature(x = "CHMfactor")}:
      see \code{\link{expand-methods}}.}
    \item{expand2}{\code{signature(x = "CHMfactor")}:
      see \code{\link{expand2-methods}}.}
    \item{image}{\code{signature(x = "CHMfactor")}:
      see \code{\link{image-methods}}.}
    \item{nnzero}{\code{signature(x = "CHMfactor")}:}
    \item{solve}{\code{signature(a = "CHMfactor", b = .)}:
      see \code{\link{solve-methods}}.}
    \item{update}{\code{signature(object = "CHMfactor")}:
      returns a copy of \code{object} with the same nonzero pattern but
      with numeric entries updated according to additional arguments
      \code{parent} and \code{mult}, where \code{parent} is (coercible to)
      a \code{\linkS4class{dsCMatrix}} or a \code{\linkS4class{dgCMatrix}}
      and \code{mult} is a numeric vector of positive length.
      
      The numeric entries are updated with those of the Cholesky factor
      of \code{F(parent)+mult[1]*I}, where the nonzero pattern
      of \code{F(parent)} \emph{must} match that of \code{S}
      if \code{object = Cholesky(S, ...)}.  Here, \code{I} is the
      identity matrix and \code{F = \link{identity}} for symmetric
      \code{parent} and \code{F = \link{tcrossprod}} otherwise.}
    \item{updown}{\code{signature(update = ., C = "mMatrix", object = "CHMfactor")}:
      see \code{\link{updown-methods}}.}
  }
}
\value{
  \code{isLDL(x)} returns \code{TRUE} or \code{FALSE}:
  \code{TRUE} if \code{x} represents a factorization \eqn{A = P' L D L' P} 
  (where \eqn{L} is unit lower triangular and \eqn{D} is diagonal),
  \code{FALSE} if \code{x} represents a factorization \eqn{A = P' L L' P}
  (where \eqn{L} is nonunit lower triangular).
}
\seealso{
  Class \code{\linkS4class{dsCMatrix}}.
  
  Generic functions \code{\link{expand2}}, \code{\link{Cholesky}},
  and \code{\link{chol}}.
}
\examples{% most are in ./Cholesky.Rd
## An example for the expand() method
n <- 1000; m <- 200; nnz <- 2000
set.seed(1)
M1 <- spMatrix(n, m,
               i = sample(n, nnz, replace = TRUE),
               j = sample(m, nnz, replace = TRUE),
               x = round(rnorm(nnz),1))
XX <- crossprod(M1) ## = M1'M1  = M M'  where M <- t(M1)
CX <- Cholesky(XX)
isLDL(CX)
str(CX) ## a "dCHMsimpl" object
r <- expand(CX)
L.P <- with(r, crossprod(L,P))  ## == L'P
PLLP <- crossprod(L.P)          ## == (L'P)' L'P == P'LL'P  = XX = M M'
b <- sample(m)
stopifnot(all.equal(PLLP, XX), % not needed: factorsCheck is off: check.attributes=FALSE),
          all(as.vector(solve(CX, b, system="P" )) == r$P \%*\% b),
          all(as.vector(solve(CX, b, system="Pt")) == t(r$P) \%*\% b) )

u1 <- update(CX, XX,    mult=pi)
u2 <- update(CX, t(M1), mult=pi) # with the original M, where XX = M M'
stopifnot(all.equal(u1,u2, tol=1e-14))

   ## [ See  help(Cholesky)  for more examples ]
   ##        -------------
}
