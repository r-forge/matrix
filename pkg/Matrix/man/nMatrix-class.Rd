\name{nMatrix-class}
\title{Virtual Class \sQuote{nMatrix} of Boolean Matrices}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{nMatrix-class}
\alias{ndenseMatrix-class}
\alias{nsparseMatrix-class}
%
\description{
  \code{nMatrix} is a virtual subclass of \code{\linkS4class{Matrix}}
  representing boolean (or \dQuote{\bold{n}onzero pattern}) matrices,
  whose entries are \code{FALSE} or \code{TRUE}.  \code{ndenseMatrix}
  and \code{nsparseMatrix} are virtual subclasses of \code{nMatrix}
  representing its intersections with \code{\linkS4class{denseMatrix}}
  and \code{\linkS4class{sparseMatrix}}.

  \code{nMatrix} objects tend to arise as the return value of
  \code{op(A)}, where \code{A} is a \code{Matrix} object and \code{op}
  is a function that preserves matrix dimensions and returns
  \code{FALSE} or \code{TRUE} by element
  (e.g., \code{op(A) := is.nan(A)}).
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}}{inherited from virtual superclass
      \code{\linkS4class{Matrix}}.}
  }
  Certain nonvirtual subclasses of \code{nMatrix} defined in package
  \pkg{Matrix} (see \sQuote{Subclasses}) also have:
  \describe{
    \item{\code{x}}{a \link{logical} vector storing matrix entries.
      Details such as \emph{which} entries are stored vary by storage
      format.  \code{NA} elements are valid and mapped to \code{TRUE},
      hence methods often must operate on not \code{x} but rather
      \code{x | is.na(x)}.}
  }
}
\section{Subclasses}{
  The nonvirtual subclasses of \code{nMatrix} are listed below with
  links to virtual superclasses defining their structure and storage
  format.  Marked (\code{*}) classes have an \code{x} slot storing
  matrix entries.  The details of the \code{x} slot vary by storage
  format.  Unmarked classes do not have an \code{x} slot.  These use
  sparse formats that store positions of nonzero entries instead of
  entries themselves.
  \tabular{rrr}{
    Nonvirtual subclass \tab Structure \tab Storage format\cr
    \code{* ngeMatrix} \tab    \code{\linkS4class{generalMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{* nsyMatrix} \tab  \code{\linkS4class{symmetricMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{* ntrMatrix} \tab \code{\linkS4class{triangularMatrix}} \tab \code{\linkS4class{unpackedMatrix}}\cr
    \code{* nspMatrix} \tab  \code{\linkS4class{symmetricMatrix}} \tab   \code{\linkS4class{packedMatrix}}\cr
    \code{* ntpMatrix} \tab \code{\linkS4class{triangularMatrix}} \tab   \code{\linkS4class{packedMatrix}}\cr
      \code{ngCMatrix} \tab    \code{\linkS4class{generalMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
      \code{nsCMatrix} \tab  \code{\linkS4class{symmetricMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
      \code{ntCMatrix} \tab \code{\linkS4class{triangularMatrix}} \tab  \code{\linkS4class{CsparseMatrix}}\cr
      \code{ngRMatrix} \tab    \code{\linkS4class{generalMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
      \code{nsRMatrix} \tab  \code{\linkS4class{symmetricMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
      \code{ntRMatrix} \tab \code{\linkS4class{triangularMatrix}} \tab  \code{\linkS4class{RsparseMatrix}}\cr
      \code{ngTMatrix} \tab    \code{\linkS4class{generalMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}\cr
      \code{nsTMatrix} \tab  \code{\linkS4class{symmetricMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}\cr
      \code{ntTMatrix} \tab \code{\linkS4class{triangularMatrix}} \tab  \code{\linkS4class{TsparseMatrix}}\cr
    \code{* ndiMatrix} \tab   \code{\linkS4class{diagonalMatrix}} \tab                \code{\var{unnamed}}\cr
      \code{indMatrix} \tab      \code{\linkS4class{indexMatrix}} \tab                \code{\var{unnamed}}
  }
  Class \code{ndiMatrix} was not defined until \pkg{Matrix} version
  1.6-2.
}
\seealso{
  \dQuote{Sister} classes \code{\linkS4class{lMatrix}},
  \code{\linkS4class{iMatrix}}, \code{\linkS4class{dMatrix}}, and
  \code{\linkS4class{zMatrix}} representing logical, integer, double,
  and complex matrices.
}
\examples{
getClass("nMatrix")

L3 <- Matrix(upper.tri(diag(3)))
L3 # an "ltCMatrix"
as(L3, "nMatrix") # -> ntC*

## similar, not using Matrix()
as(upper.tri(diag(3)), "nMatrix")# currently "ngTMatrix"

showClass("ndenseMatrix")
as(diag(3) > 0, "ndenseMatrix")# -> "nge"

(m <- Matrix(c(0,0,2:0), 3,5, dimnames=list(LETTERS[1:3],NULL)))
## ``extract the nonzero-pattern of (m) into an nMatrix'':
nm <- as(m, "nsparseMatrix") ## -> will be a "ngCMatrix"
str(nm) # no 'x' slot
nnm <- !nm # no longer sparse
## consistency check:
stopifnot(xor(as( nm, "matrix"),
              as(nnm, "matrix")))

## low-level way of adding "non-structural zeros" :
nnm <- as(nnm, "lsparseMatrix") # "lgCMatrix"
nnm@x[2:4] <- c(FALSE, NA, NA)
nnm
as(nnm, "nMatrix") # NAs *and* non-structural 0  |--->  'TRUE'

data(KNex, package = "Matrix")
nmm <- as(KNex $ mm, "nMatrix")
str(xlx <- crossprod(nmm))# "nsCMatrix"
stopifnot(isSymmetric(xlx))
image(xlx, main=paste("crossprod(nmm) : Sparse", class(xlx)))
}
