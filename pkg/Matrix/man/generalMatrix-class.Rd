\name{generalMatrix-class}
\title{Virtual Class \sQuote{generalMatrix} of General Matrices}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{generalMatrix-class}
%
\description{
  \code{generalMatrix} is a virtual subclass of
  \code{\linkS4class{Matrix}} representing general (unstructured,
  possibly nonsquare) matrices.  Nonvirtual subclasses of
  \code{generalMatrix} impose no constraints on the entries of the
  represented matrix beyond the constraint imposed by the data type.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}}{inherited from virtual superclass
      \code{\linkS4class{Matrix}}.}
    \item{\code{factors}}{a list of
      \code{\linkS4class{MatrixFactorization}} objects caching
      factorizations of the represented matrix.  It is updated
      \dQuote{automagically} by methods for generic functions that
      compute matrix factorizations, such as \code{\linkS4methods{lu}}.
      The caching mechanism is not idiomatic as it circumvents \R{}'s
      usual copy semantics; it may be deprecated in a future version of
      \pkg{Matrix}, hence code intended to be forwards compatible will
      not access the \code{factors} slot.}
  }
}
\section{Subclasses}{
  The direct, nonvirtual subclasses of \code{generalMatrix} defined in
  package \pkg{Matrix} are listed below with links to virtual
  superclasses defining their data type and storage format.
  \tabular{rrr}{
    Nonvirtual subclass \tab Data type \tab Storage format\cr
    \code{ngeMatrix} \tab \code{\linkS4class{nMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{lgeMatrix} \tab \code{\linkS4class{lMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{igeMatrix} \tab \code{\linkS4class{iMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{dgeMatrix} \tab \code{\linkS4class{dMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{zgeMatrix} \tab \code{\linkS4class{zMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{ngCMatrix} \tab \code{\linkS4class{nMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
    \code{lgCMatrix} \tab \code{\linkS4class{lMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
    \code{igCMatrix} \tab \code{\linkS4class{iMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
    \code{dgCMatrix} \tab \code{\linkS4class{dMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
    \code{zgCMatrix} \tab \code{\linkS4class{zMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
    \code{ngRMatrix} \tab \code{\linkS4class{nMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
    \code{lgRMatrix} \tab \code{\linkS4class{lMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
    \code{igRMatrix} \tab \code{\linkS4class{iMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
    \code{dgRMatrix} \tab \code{\linkS4class{dMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
    \code{zgRMatrix} \tab \code{\linkS4class{zMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
    \code{ngTMatrix} \tab \code{\linkS4class{nMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}\cr
    \code{lgTMatrix} \tab \code{\linkS4class{lMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}\cr
    \code{igTMatrix} \tab \code{\linkS4class{iMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}\cr
    \code{dgTMatrix} \tab \code{\linkS4class{dMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}\cr
    \code{zgTMatrix} \tab \code{\linkS4class{zMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}
  }
  The subclasses of \code{iMatrix} and \code{zMatrix} were not defined
  until \pkg{Matrix} version 1.8-0.
}
\seealso{
  Related classes
  \code{\linkS4class{symmetricMatrix}},
  \code{\linkS4class{posdefMatrix}},
  \code{\linkS4class{triangularMatrix}},
  \code{\linkS4class{diagonalMatrix}}, and
  \code{\linkS4class{indexMatrix}} representing
  Hermitian or symmetric, positive semidefinite, triangular, diagonal,
  and index matrices.  Generic functions \code{\linkS4methods{Schur}},
  \code{\linkS4methods{lu}}, and \code{\linkS4methods{qr}} for Schur,
  LU, and QR factorization of general matrices.
}
