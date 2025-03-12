\name{zMatrix-class}
\title{Virtual Class \sQuote{zMatrix} of Complex Matrices}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{zMatrix-class}
\alias{zdenseMatrix-class}
\alias{zsparseMatrix-class}
%
\description{
  \code{zMatrix} is a virtual subclass of \code{\linkS4class{Matrix}}
  representing complex matrices.  \code{zdenseMatrix} and
  \code{zsparseMatrix} are virtual subclasses of \code{zMatrix}
  representing its intersections with \code{\linkS4class{denseMatrix}}
  and \code{\linkS4class{sparseMatrix}}.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}}{inherited from virtual superclass
      \code{\linkS4class{Matrix}}.}
    \item{\code{x}}{a \link{complex} vector storing matrix entries.
      Details such as \emph{which} entries are stored vary by storage
      format.}
  }
}
\section{Subclasses}{
  The nonvirtual subclasses of \code{zMatrix} are listed below with
  links to virtual superclasses defining their structure and storage
  format.
  \tabular{rrr}{
    Nonvirtual subclass \tab Structure \tab Storage format\cr
    \code{zgeMatrix} \tab    \code{\linkS4class{generalMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{zsyMatrix} \tab  \code{\linkS4class{symmetricMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{zpoMatrix} \tab     \code{\linkS4class{posdefMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{ztrMatrix} \tab \code{\linkS4class{triangularMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{zspMatrix} \tab  \code{\linkS4class{symmetricMatrix}} \tab   \code{\linkS4class{packedMatrix}}\cr
    \code{zppMatrix} \tab     \code{\linkS4class{posdefMatrix}} \tab   \code{\linkS4class{packedMatrix}}\cr
    \code{ztpMatrix} \tab \code{\linkS4class{triangularMatrix}} \tab   \code{\linkS4class{packedMatrix}}\cr
    \code{zgCMatrix} \tab    \code{\linkS4class{generalMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
    \code{zsCMatrix} \tab  \code{\linkS4class{symmetricMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
    \code{zpCMatrix} \tab     \code{\linkS4class{posdefMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
    \code{ztCMatrix} \tab \code{\linkS4class{triangularMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
    \code{zgRMatrix} \tab    \code{\linkS4class{generalMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
    \code{zsRMatrix} \tab  \code{\linkS4class{symmetricMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
    \code{zpRMatrix} \tab     \code{\linkS4class{posdefMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
    \code{ztRMatrix} \tab \code{\linkS4class{triangularMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
    \code{zgTMatrix} \tab    \code{\linkS4class{generalMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}\cr
    \code{zsTMatrix} \tab  \code{\linkS4class{symmetricMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}\cr
    \code{zpTMatrix} \tab     \code{\linkS4class{posdefMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}\cr
    \code{ztTMatrix} \tab \code{\linkS4class{triangularMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}\cr
    \code{zdiMatrix} \tab   \code{\linkS4class{diagonalMatrix}} \tab                \code{\var{unnamed}}
  }
  Subclasses of \code{zMatrix} were not defined until \pkg{Matrix}
  version 1.8-0.
}
\seealso{
  \dQuote{Sister} classes \code{\linkS4class{nMatrix}},
  \code{\linkS4class{lMatrix}}, \code{\linkS4class{iMatrix}}, and
  \code{\linkS4class{dMatrix}} representing boolean, logical, integer,
  and double precision matrices.
}
\examples{
## TODO
}
