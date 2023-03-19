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
  \item{code}{a string whose first three characters specify the class
    of the result.}
  \item{kind}{a string (\code{"."}, \code{"n"}, \code{"l"}, or
    \code{"d"}) specifying the \dQuote{kind} of the result.
    \code{"."} indicates that the kind of \code{from} should
    be preserved.
    \code{"n"} indicates that the result should inherit from
    \code{"\linkS4class{nMatrix}"} (and so on).}
  \item{uplo}{a string (\code{"U"} or \code{"L"}) indicating
    whether the result should store the upper or lower triangle
    of \code{from}.  The elements of \code{from} in the opposite
    triangle are ignored.}
  \item{diag}{a string (\code{"N"} or \code{"U"}) indicating
    whether the result (if triangular) should be formally nonunit
    or unit triangular.  In the unit triangular case, the diagonal
    elements of \code{from} are ignored.}
  \item{drop0}{a logical.  If \code{TRUE}, then nonstructural
    zeros in \code{from} are dropped.}
  \item{repr}{a string (\code{"C"}, \code{"R"}, or \code{"T"})
    specifying the storage of the result as
    \code{"\linkS4class{CsparseMatrix}"},
    \code{"\linkS4class{RsparseMatrix}"}, or
    \code{"\linkS4class{TsparseMatrix}"}.}
  \item{packed}{a logical.  If \code{TRUE} and \code{from} is
    formally triangular or symmetric, then the result will have
    \dQuote{packed} storage and inherit from
    \code{"\linkS4class{packedMatrix}"} rather than
    \code{"\linkS4class{unpackedMatrix}"}.}
  \item{Csparse}{a logical.  If \code{FALSE}, then the result
    will inherit from \code{"\linkS4class{RsparseMatrix}"}
    rather than \code{"\linkS4class{CsparseMatrix}"}.}
  \item{\dots}{optional arguments passed to \code{\link{isTriangular}}
    or \code{\link{isSymmetric}}.}
}
\details{
  \code{.M2tri}, \code{.M2sym}, and \code{.M2diag} coerce their
  first argument to
  \code{"\linkS4class{triangularMatrix}"},
  \code{"\linkS4class{symmetricMatrix}"}, and
  \code{"\linkS4class{diagonalMatrix}"}, respectively,
  after checking that it has the desired structure.

  \code{.m2dense} and \code{.m2sparse} coerce their first argument,
  which must be a traditional matrix or vector,
  to \code{"\linkS4class{denseMatrix}"} and
  \code{"\link[=CsparseMatrix-class]{[CRT]sparseMatrix}"}, respectively.
  The actual class of the result is determined by \code{code},
  which must match the pattern
  \code{"^[.nld](ge|tr|sy|tp|sp)"} (for \code{.mdense}) or
  \code{"^[.nld][gts][CRT]"} (for \code{.m2sparse}),
  where \code{"."} in the first position
  means \code{"l"} for the logical type
  and \code{"d"} for the integer and double types.

  \code{.dense2m} and \code{.sparse2m} perform the reverse coercion,
  returning a traditional matrix.  \code{.diag2m} works similarly
  for \code{"\linkS4class{diagonalMatrix}"} arguments. 
}
\examples{
## TODO
}
