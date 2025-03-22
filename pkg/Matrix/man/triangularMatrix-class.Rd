\name{triangularMatrix-class}
\title{Virtual Class \sQuote{triangularMatrix} of Triangular Matrices}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{triangularMatrix-class}
%
\description{
  \code{triangularMatrix} is a virtual subclass of
  \code{\linkS4class{Matrix}} representing triangular matrices.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}}{inherited from virtual superclass
      \code{\linkS4class{Matrix}}.  \code{Dim[1]} and \code{Dim[2]} must
      be equal.}
    \item{\code{uplo}, \code{diag}}{.} % temporary
  }
  All direct, nonvirtual subclasses of \code{triangularMatrix} defined
  in package \pkg{Matrix} have two additional slots:
  \describe{
    \item{\code{uplo}}{a character string, either \code{"U"} or
      \code{"L"}.  \code{"U"} indicates that the represented matrix is
      upper triangular and that only the upper triangular entries are
      stored.  \code{"L"} indicates that the represented matrix is lower
      triangular and that only the lower triangular entries are stored.}
    \item{\code{diag}}{a character string, either \code{"N"} or
      \code{"U"}.  \code{"N"} indicates that the diagonal entries of
      the represented matrix are stored.  \code{"U"} indicates that the
      represented matrix is unit triangular and that the diagonal
      entries, which are all equal to 1, are not stored.}
  }
}
\section{Subclasses}{
  The direct, nonvirtual subclasses of \code{triangularMatrix} defined
  in package \pkg{Matrix} are listed below with links to virtual
  superclasses defining their data type and storage format.
  \tabular{rrr}{
    Nonvirtual subclass \tab Data type \tab Storage format\cr
    \code{ntrMatrix} \tab \code{\linkS4class{nMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{ltrMatrix} \tab \code{\linkS4class{lMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{itrMatrix} \tab \code{\linkS4class{iMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{dtrMatrix} \tab \code{\linkS4class{dMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{ztrMatrix} \tab \code{\linkS4class{zMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{ntpMatrix} \tab \code{\linkS4class{nMatrix}} \tab   \code{\linkS4class{packedMatrix}}\cr
    \code{ltpMatrix} \tab \code{\linkS4class{lMatrix}} \tab   \code{\linkS4class{packedMatrix}}\cr
    \code{itpMatrix} \tab \code{\linkS4class{iMatrix}} \tab   \code{\linkS4class{packedMatrix}}\cr
    \code{dtpMatrix} \tab \code{\linkS4class{dMatrix}} \tab   \code{\linkS4class{packedMatrix}}\cr
    \code{ztpMatrix} \tab \code{\linkS4class{zMatrix}} \tab   \code{\linkS4class{packedMatrix}}\cr
    \code{ntCMatrix} \tab \code{\linkS4class{nMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
    \code{ltCMatrix} \tab \code{\linkS4class{lMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
    \code{itCMatrix} \tab \code{\linkS4class{iMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
    \code{dtCMatrix} \tab \code{\linkS4class{dMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
    \code{ztCMatrix} \tab \code{\linkS4class{zMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
    \code{ntRMatrix} \tab \code{\linkS4class{nMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
    \code{ltRMatrix} \tab \code{\linkS4class{lMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
    \code{itRMatrix} \tab \code{\linkS4class{iMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
    \code{dtRMatrix} \tab \code{\linkS4class{dMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
    \code{ztRMatrix} \tab \code{\linkS4class{zMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
    \code{ntTMatrix} \tab \code{\linkS4class{nMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}\cr
    \code{ltTMatrix} \tab \code{\linkS4class{lMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}\cr
    \code{itTMatrix} \tab \code{\linkS4class{iMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}\cr
    \code{dtTMatrix} \tab \code{\linkS4class{dMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}\cr
    \code{ztTMatrix} \tab \code{\linkS4class{zMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}
  }
  The subclasses of \code{iMatrix} and \code{zMatrix} were not defined
  until \pkg{Matrix} version 1.8-0.
}
\seealso{
  Related classes
  \code{\linkS4class{generalMatrix}},
  \code{\linkS4class{symmetricMatrix}},
  \code{\linkS4class{posdefMatrix}},
  \code{\linkS4class{diagonalMatrix}}, and
  \code{\linkS4class{indexMatrix}} representing
  general (unstructured, possibly nonsquare), Hermitian or symmetric,
  positive semidefinite, diagonal, and index matrices.  Generic function
  \code{\linkS4methods{isTriangular}} for testing for testing if a
  matrix is upper or lower triangular.  Generic function
  \code{\linkS4methods{forceTriangular}} for inducing upper or lower
  triangular structure.  Generic functions \code{\linkS4methods{band}},
  \code{\linkS4methods{triu}}, and \code{\linkS4methods{tril}} for
  inducing banded, upper trapezoidal, and lower trapezoidal structure.
}
\examples{
showClass("triangularMatrix")

## The names of direct subclasses:
scl <- getClass("triangularMatrix")@subclasses
directly <- sapply(lapply(scl, slot, "by"), length) == 0
names(scl)[directly]

(m <- matrix(c(5, 1, 0, 3), 2))
as(m, "triangularMatrix")
}
