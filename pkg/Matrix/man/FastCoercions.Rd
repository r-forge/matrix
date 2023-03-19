\name{FastCoercions}
\title{Low Level Coercion Functions}
%
\alias{FastCoercions}
%
\alias{.CR2RC}
\alias{.CR2T}
\alias{.M2diag}
\alias{.M2sym}
\alias{.M2tri}
\alias{.T2CR}
\alias{.dense2g}
\alias{.dense2kind}
\alias{.dense2m}
\alias{.dense2sparse}
\alias{.dense2v}
\alias{.diag2dense}
\alias{.diag2sparse}
\alias{.m2dense}
\alias{.m2sparse}
\alias{.sparse2dense}
\alias{.sparse2g}
\alias{.sparse2kind}
\alias{.sparse2m}
\alias{.sparse2v}
\alias{.tCR2RC}
%
\description{
  \dQuote{Semi-API} functions used internally by \pkg{Matrix},
  often to bypass S4 dispatch and avoid the associated overhead.
  These are exported to provide this capability to expert users.
  Typical users should continue to rely on S4 generic functions
  to dispatch suitable methods.
}
\usage{
.M2tri(from, ...)
.M2sym(from, ...)
.M2diag(from)

.m2dense(from, code, uplo = "U", diag = "N")
.m2sparse(from, code, uplo = "U", diag = "N")

.dense2m(from)
.sparse2m(from)

.dense2v(from)
.sparse2v(from)

.dense2kind(from, kind)
.sparse2kind(from, kind, drop0 = FALSE)

.dense2g(from, kind = ".")
.sparse2g(from)

.dense2sparse(from, repr = "C")
.sparse2dense(from, packed = FALSE)

.diag2dense(from, code, uplo = "U")
.diag2sparse(from, code, uplo = "U", drop0 = TRUE)

.CR2T(from)
.T2CR(from, Csparse = TRUE)

.CR2RC(from)
.tCR2RC(from)
}
\arguments{
  \item{from}{a \code{\linkS4class{Matrix}}, matrix, or vector.}
  \item{code}{a string whose first three characters specify the class of
    the result.  It should match the pattern
    \code{"^[.nld](ge|tr|sy|tp|sp)"} for \code{.*2dense} and
    \code{"^[.nld][gts][CRT]"} for \code{.*2sparse},
    where \code{"."} in the first position is equivalent to \code{"l"}
    for logical arguments and \code{"d"} for numeric arguments.}
  \item{kind}{a string (\code{"."}, \code{"n"}, \code{"l"}, or
    \code{"d"}) specifying the \dQuote{kind} of the result.  \code{"."}
    indicates that the kind of \code{from} should be preserved.
    \code{"n"} indicates that the result should inherit from
    \code{\linkS4class{nMatrix}} (and so on).}
  \item{uplo}{a string (\code{"U"} or \code{"L"}) indicating whether the
    result should store the upper or lower triangle of \code{from}.  The
    elements of \code{from} in the opposite triangle are ignored.}
  \item{diag}{a string (\code{"N"} or \code{"U"}) indicating whether the
    result (if triangular) should be formally nonunit or unit
    triangular.  In the unit triangular case, the diagonal elements of
    \code{from} are ignored.}
  \item{drop0}{a logical.  If \code{TRUE}, then nonstructural zeros in
    \code{from} are dropped.}
  \item{repr}{a string (\code{"C"}, \code{"R"}, or \code{"T"})
    specifying the storage of the result as
    \code{\linkS4class{CsparseMatrix}},
    \code{\linkS4class{RsparseMatrix}}, or
    \code{\linkS4class{TsparseMatrix}}.}
  \item{packed}{a logical.  If \code{TRUE} and \code{from} is formally
    triangular or symmetric, then the result will have \dQuote{packed}
    storage and inherit from \code{\linkS4class{packedMatrix}} rather
    than \code{\linkS4class{unpackedMatrix}}.}
  \item{Csparse}{a logical.  If \code{FALSE}, then the result will
    inherit from \code{\linkS4class{RsparseMatrix}} rather than
    \code{\linkS4class{CsparseMatrix}}.}
  \item{\dots}{optional arguments passed to \code{\link{isTriangular}}
    or \code{\link{isSymmetric}}.}
}
\details{
  Functions with names of the form \code{.<A>2<B>} implement coercions
  from virtual class A to the \dQuote{nearest} non-virtual subclass of
  virtual class B, where the virtual classes are abbreviated as follows:

  \describe{
    \item{M}{\code{\linkS4class{Matrix}}, matrix, or vector}
    \item{m}{matrix}
    \item{v}{vector}
    \item{g}{\code{\linkS4class{generalMatrix}}}
    \item{C}{\code{\linkS4class{CsparseMatrix}}}
    \item{R}{\code{\linkS4class{RsparseMatrix}}}
    \item{T}{\code{\linkS4class{TsparseMatrix}}}
    \item{dense}{\code{\linkS4class{denseMatrix}}}
    \item{sparse}{\code{\linkS4class{CsparseMatrix}},
      \code{\linkS4class{RsparseMatrix}}, or
      \code{\linkS4class{TsparseMatrix}}}
    \item{tri}{\code{\linkS4class{triangularMatrix}}}
    \item{sym}{\code{\linkS4class{symmetricMatrix}}}
    \item{diag}{\code{\linkS4class{diagonalMatrix}}}
  }
  
  The names should be seen as guides, rather than as an exact
  description of behaviour.  For example, \code{.m2dense} and
  \code{.m2sparse} accept vectors in addition to matrices.
  
  \code{.CR2T} and \code{.T2CR} coerce between
  \code{\linkS4class{TsparseMatrix}} and the union of
  \code{\linkS4class{CsparseMatrix}} and
  \code{\linkS4class{RsparseMatrix}}.

  \code{.CR2RC} and \code{.tCR2RC} coerce between
  \code{\linkS4class{CsparseMatrix}} and
  \code{\linkS4class{RsparseMatrix}}.  Conceptually, the latter
  performs the coercion on the transpose of its argument.
  That is, \code{.tCR2RC(from)} is equivalent to but much more
  efficient than \code{.CR2RC(t(from))} and \code{t(.CR2RC(from))}.
  
  \code{.M2tri}, \code{.M2sym}, and \code{.M2diag} can be seen as
  drop-in replacements for \code{as(., "*Matrix")}, but allowing users
  to pass optional arguments to the structure-checking functions.
}
\examples{
## TODO
}
