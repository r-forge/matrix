\name{symmetricMatrix-class}
\title{Virtual Class \sQuote{symmetricMatrix} of Hermitian or Symmetric Matrices}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{symmetricMatrix-class}
%
\description{
  \code{symmetricMatrix} is a virtual subclass of
  \code{\linkS4class{Matrix}} representing Hermitian or symmetric
  matrices.  Objects inheriting from \code{symmetricMatrix} have
  symmetric dimension names.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}}{inherited from virtual superclass
      \code{\linkS4class{Matrix}}.  \code{Dim[1]} and \code{Dim[2]} must
      be equal.}
    \item{\code{Dimnames}}{inherited from virtual superclass
      \code{\linkS4class{Matrix}}.  The actual dimension names of
      objects \code{x} are obtained by \code{dimnames(x)}, not by
      \code{x@Dimnames}.  The method for \code{dimnames} symmetrizes the
      \code{Dimnames} slot, which can be asymmetric, and returns the
      result; see \code{\link{dimnames-methods}}.  The validity method
      for class \code{symmetricMatrix} may be changed in a future
      version of \pkg{Matrix} to require a symmetric \code{Dimnames}
      slot, so that the actual dimension names of objects \code{x} can
      be obtained directly as \code{x@Dimnames}.}
    \item{\code{factors}}{a list of
      \code{\linkS4class{MatrixFactorization}} objects caching
      factorizations of the represented matrix.  It is updated
      \dQuote{automagically} by methods for generic functions that
      compute matrix factorizations, such as \code{\linkS4methods{lu}}.
      The caching mechanism is not idiomatic as it circumvents \R{}'s
      usual copy semantics; it may be deprecated in a future version of
      \pkg{Matrix}, hence code intended to be forwards compatible will
      not access the \code{factors} slot.}
    \item{\code{uplo}}{.} % temporary
  }
  All direct, nonvirtual subclasses of \code{symmetricMatrix} defined in
  package \pkg{Matrix} have an additional slot:
  \describe{
    \item{\code{uplo}}{a character string, either \code{"U"} or
      \code{"L"}.  \code{"U"} indicates that only the upper triangular
      entries of the represented matrix are stored.  \code{"L"}
      indicates that only the lower triangular entries of the
      represented matrix are stored.}
  }
}
\section{Subclasses}{
  The direct, nonvirtual subclasses of \code{symmetricMatrix} defined in
  package \pkg{Matrix} are listed below with links to virtual
  superclasses defining their data type and storage format.
  \tabular{rrr}{
    Nonvirtual subclass \tab Data type \tab Storage format\cr
    \code{nsyMatrix} \tab \code{\linkS4class{nMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{lsyMatrix} \tab \code{\linkS4class{lMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{isyMatrix} \tab \code{\linkS4class{iMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{dsyMatrix} \tab \code{\linkS4class{dMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{zsyMatrix} \tab \code{\linkS4class{zMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{nspMatrix} \tab \code{\linkS4class{nMatrix}} \tab   \code{\linkS4class{packedMatrix}}\cr
    \code{lspMatrix} \tab \code{\linkS4class{lMatrix}} \tab   \code{\linkS4class{packedMatrix}}\cr
    \code{ispMatrix} \tab \code{\linkS4class{iMatrix}} \tab   \code{\linkS4class{packedMatrix}}\cr
    \code{dspMatrix} \tab \code{\linkS4class{dMatrix}} \tab   \code{\linkS4class{packedMatrix}}\cr
    \code{zspMatrix} \tab \code{\linkS4class{zMatrix}} \tab   \code{\linkS4class{packedMatrix}}\cr
    \code{ngCMatrix} \tab \code{\linkS4class{nMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
    \code{lsCMatrix} \tab \code{\linkS4class{lMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
    \code{isCMatrix} \tab \code{\linkS4class{iMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
    \code{dsCMatrix} \tab \code{\linkS4class{dMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
    \code{zsCMatrix} \tab \code{\linkS4class{zMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
    \code{nsRMatrix} \tab \code{\linkS4class{nMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
    \code{lsRMatrix} \tab \code{\linkS4class{lMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
    \code{isRMatrix} \tab \code{\linkS4class{iMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
    \code{dsRMatrix} \tab \code{\linkS4class{dMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
    \code{zsRMatrix} \tab \code{\linkS4class{zMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
    \code{nsTMatrix} \tab \code{\linkS4class{nMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}\cr
    \code{lsTMatrix} \tab \code{\linkS4class{lMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}\cr
    \code{isTMatrix} \tab \code{\linkS4class{iMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}\cr
    \code{dsTMatrix} \tab \code{\linkS4class{dMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}\cr
    \code{zsTMatrix} \tab \code{\linkS4class{zMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}
  }
  The subclasses of \code{iMatrix} and \code{zMatrix} were not defined
  until \pkg{Matrix} version 1.8-0.
}
\seealso{
  Related classes
  \code{\linkS4class{generalMatrix}},
  \code{\linkS4class{posdefMatrix}},
  \code{\linkS4class{triangularMatrix}},
  \code{\linkS4class{diagonalMatrix}}, and
  \code{\linkS4class{indexMatrix}} representing
  general (unstructured, possibly nonsquare), positive semidefinite,
  triangular, diagonal, and index matrices.  Generic function
  \code{\linkS4methods{isSymmetric}} for testing for testing if a matrix
  is (possibly \emph{approximately}) Hermitian or symmetric.  Generic
  function \code{\linkS4methods{forceSymmetric}} for inducing Hermitian
  or symmetric structure.  Generic function
  \code{\linkS4methods{symmpart}} for getting the Hermitian or symmetric
  part of a square matrix.  Generic function
  \code{\linkS4methods{skewpart}} for getting the skew-Hermitian or
  skew-symmetric part of a square matrix.  Generic function
  \code{\linkS4methods{BunchKaufman}} for Bunch-Kaufman factorization of
  Hermitian or symmetric matrices.
}
\examples{
## An example about the symmetric Dimnames:
sy <- sparseMatrix(i= c(2,4,3:5), j= c(4,7:5,5), x = 1:5, dims = c(7,7),
                   symmetric=TRUE, dimnames = list(NULL, letters[1:7]))
sy # shows symmetrical dimnames
sy@Dimnames  # internally only one part is stored
dimnames(sy) # both parts - as sy *is* symmetrical
\dontshow{
local({ nm <- letters[1:7]
  stopifnot(identical(dimnames(sy), list(  nm, nm)),
	    identical(sy@Dimnames , list(NULL, nm)))
})
}%dont
showClass("symmetricMatrix")

## The names of direct subclasses:
scl <- getClass("symmetricMatrix")@subclasses
directly <- sapply(lapply(scl, slot, "by"), length) == 0
names(scl)[directly]

## Methods -- applicaple to all subclasses above:
showMethods(classes = "symmetricMatrix")
}
