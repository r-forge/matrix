\name{lsparseMatrix-class}
\title{Sparse logical matrices}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{lsparseMatrix-class}
\alias{lsparseMatrix-classes} % used by package 'enhancerHomologSearch'
%
\alias{coerce,matrix,lsparseMatrix-method}
\alias{coerce,vector,lsparseMatrix-method}
%
\description{The \code{lsparseMatrix} class is a virtual class
  of logical sparse matrices, i.e., sparse matrices with entries
  \code{TRUE}, \code{FALSE}, or \code{NA}.

  These can be stored in the \dQuote{triplet} form (class
  \code{\linkS4class{TsparseMatrix}}, subclasses \code{lgTMatrix},
  \code{lsTMatrix}, and \code{ltTMatrix}) or in compressed
  column-oriented form (class \code{\linkS4class{CsparseMatrix}},
  subclasses \code{lgCMatrix}, \code{lsCMatrix}, and \code{ltCMatrix})
  or--\emph{rarely}--in compressed row-oriented form (class
  \code{\linkS4class{RsparseMatrix}}, subclasses \code{lgRMatrix},
  \code{lsRMatrix}, and \code{ltRMatrix}).  The second letter in the
  name of these non-virtual classes indicates \code{g}eneral,
  \code{s}ymmetric, or \code{t}riangular.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("lgCMatrix",
    ...)} and so on.  More frequently objects are created by coercion of
  a numeric sparse matrix to the logical form, e.g. in an expression
  \code{x != 0}.

  The logical form is also used in the symbolic analysis phase
  of an algorithm involving sparse matrices.  Such algorithms often
  involve two phases: a symbolic phase wherein the positions of the
  non-zeros in the result are determined and a numeric phase wherein the
  actual results are calculated.  During the symbolic phase only the
  positions of the non-zero elements in any operands are of interest,
  hence any numeric sparse matrices can be treated as logical sparse
  matrices.
}
\details{
  Note that triplet stored (\code{\linkS4class{TsparseMatrix}}) matrices
  such as \code{lgTMatrix} may contain duplicated pairs of indices
  \eqn{(i,j)} as for the corresponding numeric class
  \code{\linkS4class{dgTMatrix}} where for such pairs, the corresponding
  \code{x} slot entries are added.  For logical matrices, the \code{x}
  entries corresponding to duplicated index pairs \eqn{(i,j)} are
  \dQuote{added} as well if the addition is defined as logical \eqn{or},
  i.e., \dQuote{\code{TRUE + TRUE  |->  TRUE}} and
  \dQuote{\code{TRUE + FALSE  |->  TRUE}}.
  Note the use of \code{\link{asUniqueT}()} for getting an internally
  unique representation without duplicated \eqn{(i,j)} entries.
}
\section{Slots}{
  \describe{
    \item{\code{x}:}{Object of class \code{"logical"}, i.e., either
      \code{TRUE}, \code{\link{NA}}, or \code{FALSE}.}
    \item{\code{Dim}:}{Object of class \code{"integer"} - the dimensions
      of the matrix.}
  }
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "dgCMatrix", to = "lgCMatrix")}}
    \item{t}{\code{signature(x = "lgCMatrix")}: returns the transpose
      of \code{x}}

    \item{which}{\code{signature(x = "lsparseMatrix")}, semantically
      equivalent to \pkg{base} function \code{\link{which}(x, arr.ind)};
      for details, see the \code{\linkS4class{lMatrix}} class documentation.}
  }
}
\seealso{
  the class \code{\linkS4class{dgCMatrix}} and \code{\linkS4class{dgTMatrix}}
}
\examples{
\dontshow{ % for R_DEFAULT_PACKAGES=NULL
library(utils, pos = "package:base", verbose = FALSE)
}
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
