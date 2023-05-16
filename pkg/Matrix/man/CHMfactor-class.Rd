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
\alias{coerce,CHMfactor,dtCMatrix-method}
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
  factorizations of \eqn{n \times n}{n-by-n} real, symmetric,
  positive semidefinite matrices, having the general form
  \deqn{P_{1} A P_{1}' = L_{1} D L_{1}' = L L'}{P1 * A * P1' = L1 * D * L1' = L * L'}
  or (equivalently)
  \deqn{A = P_{1}' L_{1} D L_{1}' P_{1} = P_{1}' L L' P_{1}}{A = P1' L1 * D * L1' * P1 = P1' * L * L' * P1}
  where
  \eqn{P_{1}}{P1} is a permutation matrix,
  \eqn{L_{1}} is a unit lower triangular matrix,
  \code{D} is a non-negative diagonal matrix, and
  \eqn{L = L_{1} \sqrt{D}}{L = L1 * sqrt(D)}.

  The implementation of class \code{CHMfactor} is based on
  CHOLMOD's C-level \code{cholmod_factor_struct}. Virtual
  subclasses \code{CHMsimpl} and \code{CHMsuper} separate
  the simplicial and supernodal variants.  These have nonvirtual
  subclasses \code{[dn]CHMsimpl} and \code{[dn]CHMsuper}, where
  prefix \samp{d} and prefix \samp{n} are reserved for numeric
  and symbolic factorizations, respectively.
}
\usage{
isLDL(x)
}
\arguments{
  \item{x}{an object inheriting from virtual class \code{CHMfactor},
    almost always the result of a call to generic function
    \code{\link{Cholesky}}.}
}
\value{
  \code{isLDL(x)} returns \code{TRUE} or \code{FALSE}:
  \code{TRUE} if \code{x} stores the diagonal entries of \eqn{D}
  and strictly lower triangular entries of \eqn{L_{1}}{L1},
  \code{FALSE} if \code{x} stores the lower triangular entries
  of \eqn{L = L_{1} \sqrt{D}}{L = L1 * sqrt(D)}.
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
    \item{\code{perm}}{a 0-based integer vector of length \code{Dim[1]}
      specifying the permutation applied to the rows and columns of the
      factorized matrix.
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
    %% MJ: Does CHOLMOD use p[length(p)] ??
    \item{\code{i}}{an integer vector of length greater than or equal
      to \code{sum(nz)} containing the row indices of nonzero entries
      in the lower triangular Cholesky factor.  These are grouped by
      column and sorted within columns, but the columns themselves
      need not be ordered monotonically.  Columns may be overallocated,
      i.e., the number of elements of \code{i} reserved for column
      \code{j} may exceed \code{nz[j]}.}
    \item{\code{prv}, \code{nxt}}{integer vectors of length
      \code{Dim[1]+2} indicating the order in which the columns of
      the lower triangular Cholesky factor are stored in \code{i}
      and \code{x}.}
    %% MJ: _if_ the columns are ordered monotonically ...
    %% For \code{j} in \code{seq_len(n <- Dim[1])}:
    %% \code{prv[1] = n+1}, \code{prv[n+2] = -1}, \code{prv[j+1] = j-1},
    %% \code{nxt[n+1] = -1}, \code{nxt[n+2] = 0}, \code{nxt[j] = j}.
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
    \item{\code{chol2inv}}{\code{signature(x = "CHMfactor")}:
      see \code{\link{chol2inv-methods}}.}
    \item{\code{coerce}}{\code{signature(from = "CHMfactor", to = "dtCMatrix")}:
      returns a \code{\linkS4class{dtCMatrix}} representing
      the lower triangular Cholesky factor \eqn{L}.}
    \item{\code{coerce}}{\code{signature(from = "CHMfactor", to = "pMatrix")}:
      returns a \code{\linkS4class{pMatrix}} representing
      the permutation matrix \eqn{P_{1}}.}
    \item{\code{determinant}}{\code{signature(from = "CHMfactor", logarithm = "logical")}:
      computes the determinant of the factorized matrix \eqn{A}
      or its logarithm.}
    \item{\code{expand}}{\code{signature(x = "CHMfactor")}:
      see \code{\link{expand-methods}}.}
    \item{\code{expand2}}{\code{signature(x = "CHMfactor")}:
      see \code{\link{expand2-methods}}.}
    \item{\code{image}}{\code{signature(x = "CHMfactor")}:
      see \code{\link{image-methods}}.}
    \item{\code{nnzero}}{\code{signature(x = "CHMfactor")}:}
    \item{\code{solve}}{\code{signature(a = "CHMfactor", b = .)}:
      see \code{\link{solve-methods}}.}
    \item{\code{update}}{\code{signature(object = "CHMfactor")}:
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
    \item{\code{updown}}{\code{signature(update = ., C = "mMatrix", object = "CHMfactor")}:
      see \code{\link{updown-methods}}.}
  }
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
