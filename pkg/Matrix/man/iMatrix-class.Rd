\name{iMatrix-class}
\title{Virtual Class \sQuote{iMatrix} of Integer Matrices}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{iMatrix-class}
\alias{idenseMatrix-class}
\alias{isparseMatrix-class}
%
\description{
  \code{iMatrix} is a virtual subclass of \code{\linkS4class{Matrix}}
  representing integer matrices, whose data are stored using numeric
  vectors of type \code{"integer"}.  \code{idenseMatrix} and
  \code{isparseMatrix} are virtual subclasses of \code{iMatrix}
  representing its intersections with \code{\linkS4class{denseMatrix}}
  and \code{\linkS4class{sparseMatrix}}.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}}{inherited from virtual superclass
      \code{\linkS4class{Matrix}}.}
    \item{\code{x}}{an \link{integer} vector storing matrix entries.
      Details such as \emph{which} entries are stored vary by storage
      format.}
  }
}
\section{Subclasses}{
  The nonvirtual subclasses of \code{iMatrix} are listed below with
  links to virtual superclasses defining their structure and storage
  format.
  \tabular{rrr}{
    Nonvirtual subclass \tab Structure \tab Storage format\cr
    \code{igeMatrix} \tab    \code{\linkS4class{generalMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{isyMatrix} \tab  \code{\linkS4class{symmetricMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{itrMatrix} \tab \code{\linkS4class{triangularMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{ispMatrix} \tab  \code{\linkS4class{symmetricMatrix}} \tab   \code{\linkS4class{packedMatrix}}\cr
    \code{itpMatrix} \tab \code{\linkS4class{triangularMatrix}} \tab   \code{\linkS4class{packedMatrix}}\cr
    \code{igCMatrix} \tab    \code{\linkS4class{generalMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
    \code{isCMatrix} \tab  \code{\linkS4class{symmetricMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
    \code{itCMatrix} \tab \code{\linkS4class{triangularMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
    \code{igRMatrix} \tab    \code{\linkS4class{generalMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
    \code{isRMatrix} \tab  \code{\linkS4class{symmetricMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
    \code{itRMatrix} \tab \code{\linkS4class{triangularMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
    \code{igTMatrix} \tab    \code{\linkS4class{generalMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}\cr
    \code{isTMatrix} \tab  \code{\linkS4class{symmetricMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}\cr
    \code{itTMatrix} \tab \code{\linkS4class{triangularMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}\cr
    \code{idiMatrix} \tab   \code{\linkS4class{diagonalMatrix}} \tab                \code{\var{unnamed}}
  }
  Subclasses of \code{iMatrix} were not defined until \pkg{Matrix}
  version 1.8-0.
}
\seealso{
  \dQuote{Sister} classes \code{\linkS4class{nMatrix}},
  \code{\linkS4class{lMatrix}}, \code{\linkS4class{dMatrix}}, and
  \code{\linkS4class{zMatrix}} representing boolean, logical, double,
  and complex matrices.
}
\examples{
## TODO
}
