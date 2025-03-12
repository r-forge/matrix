\name{dMatrix-class}
\title{Virtual Class \sQuote{dMatrix} of Double Precision Matrices}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{dMatrix-class}
\alias{ddenseMatrix-class}
\alias{dsparseMatrix-class}
%
\description{
  \code{dMatrix} is a virtual subclass of \code{\linkS4class{Matrix}}
  representing double precision matrices, whose data are stored using
  numeric vectors of type \code{"double"}.  \code{ddenseMatrix} and
  \code{dsparseMatrix} are virtual subclasses of \code{dMatrix}
  representing its intersections with \code{\linkS4class{denseMatrix}}
  and \code{\linkS4class{sparseMatrix}}.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}}{inherited from virtual superclass
      \code{\linkS4class{Matrix}}.}
    \item{\code{x}}{a \link{double} vector storing matrix entries.
      Details such as \emph{which} entries are stored vary by storage
      format.}
  }
}
\section{Subclasses}{
  The nonvirtual subclasses of \code{dMatrix} are listed below with
  links to virtual superclasses defining their structure and storage
  format.
  \tabular{rrr}{
    Nonvirtual subclass \tab Structure \tab Storage format\cr
    \code{dgeMatrix} \tab    \code{\linkS4class{generalMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{dsyMatrix} \tab  \code{\linkS4class{symmetricMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{dpoMatrix} \tab     \code{\linkS4class{posdefMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{dtrMatrix} \tab \code{\linkS4class{triangularMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{dspMatrix} \tab  \code{\linkS4class{symmetricMatrix}} \tab   \code{\linkS4class{packedMatrix}}\cr
    \code{dppMatrix} \tab     \code{\linkS4class{posdefMatrix}} \tab   \code{\linkS4class{packedMatrix}}\cr
    \code{dtpMatrix} \tab \code{\linkS4class{triangularMatrix}} \tab   \code{\linkS4class{packedMatrix}}\cr
    \code{dgCMatrix} \tab    \code{\linkS4class{generalMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
    \code{dsCMatrix} \tab  \code{\linkS4class{symmetricMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
    \code{dpCMatrix} \tab     \code{\linkS4class{posdefMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
    \code{dtCMatrix} \tab \code{\linkS4class{triangularMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
    \code{dgRMatrix} \tab    \code{\linkS4class{generalMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
    \code{dsRMatrix} \tab  \code{\linkS4class{symmetricMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
    \code{dpRMatrix} \tab     \code{\linkS4class{posdefMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
    \code{dtRMatrix} \tab \code{\linkS4class{triangularMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
    \code{dgTMatrix} \tab    \code{\linkS4class{generalMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}\cr
    \code{dsTMatrix} \tab  \code{\linkS4class{symmetricMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}\cr
    \code{dpTMatrix} \tab     \code{\linkS4class{posdefMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}\cr
    \code{dtTMatrix} \tab \code{\linkS4class{triangularMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}\cr
    \code{ddiMatrix} \tab   \code{\linkS4class{diagonalMatrix}} \tab                \code{\var{unnamed}}
  }
  Classes \code{dpCMatrix}, \code{dpRMatrix}, and \code{dpTMatrix}
  were not defined until \pkg{Matrix} version 1.8-0.
}
\seealso{
  \dQuote{Sister} classes \code{\linkS4class{nMatrix}},
  \code{\linkS4class{lMatrix}}, \code{\linkS4class{iMatrix}}, and
  \code{\linkS4class{zMatrix}} representing boolean, logical, integer,
  and complex matrices.
}
\examples{
showClass("dMatrix")

set.seed(101)
round(Matrix(rnorm(28), 4, 7), 2)
M <- Matrix(rlnorm(56, sd = 10), 4, 14)
(M. <- zapsmall(M))
table(as.logical(M. == 0))

showClass("ddenseMatrix")
showMethods(class = "ddenseMatrix", where = "package:Matrix")

showClass("dsparseMatrix")
}
