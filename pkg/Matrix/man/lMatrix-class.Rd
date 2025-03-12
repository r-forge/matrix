\name{lMatrix-class}
\title{Virtual Class \sQuote{lMatrix} of Logical Matrices}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{lMatrix-class}
\alias{ldenseMatrix-class}
\alias{lsparseMatrix-class}
%
\description{
  \code{lMatrix} is a virtual subclass of \code{\linkS4class{Matrix}}
  representing logical matrices, whose entries are
  \code{FALSE}, \code{TRUE}, or \code{NA}.  \code{ldenseMatrix} and
  \code{lsparseMatrix} are virtual subclasses of \code{lMatrix}
  representing its intersections with \code{\linkS4class{denseMatrix}}
  and \code{\linkS4class{sparseMatrix}}.

  \code{lMatrix} objects tend to arise as the return value of
  \code{op(A)}, where \code{A} is a \code{Matrix} object and \code{op}
  is a function that preserves matrix dimensions and returns
  \code{FALSE}, \code{TRUE}, or \code{NA} by element
  (e.g., \code{op(A) := A == 0} or \code{!A} or \code{A | B} for some
  conformable \code{B}).
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}}{inherited from virtual superclass
      \code{\linkS4class{Matrix}}.}
    \item{\code{x}}{a \link{logical} vector storing matrix entries.
      Details such as \emph{which} entries are stored vary by storage
      format.}
  }
}
\section{Subclasses}{
  The nonvirtual subclasses of \code{lMatrix} are listed below with
  links to virtual superclasses defining their structure and storage
  format.
  \tabular{rrr}{
    Nonvirtual subclass \tab Structure \tab Storage format\cr
    \code{lgeMatrix} \tab    \code{\linkS4class{generalMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{lsyMatrix} \tab  \code{\linkS4class{symmetricMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{ltrMatrix} \tab \code{\linkS4class{triangularMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{lspMatrix} \tab  \code{\linkS4class{symmetricMatrix}} \tab   \code{\linkS4class{packedMatrix}}\cr
    \code{ltpMatrix} \tab \code{\linkS4class{triangularMatrix}} \tab   \code{\linkS4class{packedMatrix}}\cr
    \code{lgCMatrix} \tab    \code{\linkS4class{generalMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
    \code{lsCMatrix} \tab  \code{\linkS4class{symmetricMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
    \code{ltCMatrix} \tab \code{\linkS4class{triangularMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
    \code{lgRMatrix} \tab    \code{\linkS4class{generalMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
    \code{lsRMatrix} \tab  \code{\linkS4class{symmetricMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
    \code{ltRMatrix} \tab \code{\linkS4class{triangularMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
    \code{lgTMatrix} \tab    \code{\linkS4class{generalMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}\cr
    \code{lsTMatrix} \tab  \code{\linkS4class{symmetricMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}\cr
    \code{ltTMatrix} \tab \code{\linkS4class{triangularMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}\cr
    \code{ldiMatrix} \tab   \code{\linkS4class{diagonalMatrix}} \tab                \code{\var{unnamed}}
  }
}
\seealso{
  \dQuote{Sister} classes \code{\linkS4class{nMatrix}},
  \code{\linkS4class{iMatrix}}, \code{\linkS4class{dMatrix}}, and
  \code{\linkS4class{zMatrix}} representing boolean, integer, double,
  and complex matrices.
}
\examples{
showClass("ldenseMatrix")

as(diag(3) > 0, "ldenseMatrix")

(m <- Matrix(c(0,0,2:0), 3,5, dimnames=list(LETTERS[1:3],NULL)))
(lm <- (m > 1)) # lgC
!lm     # no longer sparse
stopifnot(is(lm,"lsparseMatrix"),
          identical(!lm, m <= 1))

data(KNex, package = "Matrix")
str(mmG.1 <- (KNex $ mm) > 0.1)# "lgC..."
table(mmG.1@x)# however with many ``non-structural zeros''
## from logical to nz_pattern -- okay when there are no NA's :
nmG.1 <- as(mmG.1, "nMatrix") # <<< has "TRUE" also where mmG.1 had FALSE
## from logical to "double"
dmG.1 <- as(mmG.1, "dMatrix") # has '0' and back:
lmG.1 <- as(dmG.1, "lMatrix")
stopifnot(identical(nmG.1, as((KNex $ mm) != 0,"nMatrix")),
          validObject(lmG.1),
          identical(lmG.1, mmG.1))

class(xnx <- crossprod(nmG.1))# "nsC.."
class(xlx <- crossprod(mmG.1))# "dsC.." : numeric
is0 <- (xlx == 0)
mean(as.vector(is0))# 99.3\% zeros: quite sparse, but
table(xlx@x == 0)# more than half of the entries are (non-structural!) 0
stopifnot(isSymmetric(xlx), isSymmetric(xnx),
          ## compare xnx and xlx : have the *same* non-structural 0s :
          sapply(slotNames(xnx),
                 function(n) identical(slot(xnx, n), slot(xlx, n))))
}
