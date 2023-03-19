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
  
  Abbreviations should be seen as guides, rather than as an exact
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
D. <- diag(x = c(1, 1, 2, 3, 5, 8))
D.0 <- Diagonal(x = c(0, 0, 0, 3, 5, 8))
S. <- toeplitz(as.double(1:6))
C. <- new("dgCMatrix", Dim = c(3L, 4L),
          p = c(0L, 1L, 1L, 1L, 3L), i = c(1L, 0L, 2L), x = c(-8, 2, 3))

stopifnot(identical(.M2tri( D.), as(D., "triangularMatrix")),
          identical(.M2sym( D.), as(D., "symmetricMatrix")),
          identical(.M2diag(D.), as(D., "diagonalMatrix")),
          identical(.sparse2kind(C., "l"),
                    as(C., "lMatrix")),
          identical(.dense2kind(.sparse2dense(C.), "l"),
                    as(as(C., "denseMatrix"), "lMatrix")),
          identical(.diag2sparse(D.0, "ntC"),
                    .dense2sparse(.diag2dense(D.0, "ntp"), "C")),
          identical(.dense2g(.diag2dense(D.0, "dsy")),
                    .sparse2dense(.sparse2g(.diag2sparse(D.0, "dsT")))),
          identical(S.,
                    .sparse2m(.m2sparse(S., ".sR"))),
          identical(S. * lower.tri(S.) + diag(1, 6L),
                    .dense2m(.m2dense(S., ".tr", "L", "U"))),
          identical(.CR2RC(C.), .T2CR(.CR2T(C.), FALSE)),
          identical(.tCR2RC(C.), .CR2RC(t(C.))))
}
