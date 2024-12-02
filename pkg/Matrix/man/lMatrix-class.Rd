\name{lMatrix-class}
\title{Virtual Class "lMatrix" of Logical Matrices}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{lMatrix-class}
\alias{ldenseMatrix-class}
\alias{lsparseMatrix-class}
%
\alias{coerce,matrix,lMatrix-method}
\alias{coerce,vector,lMatrix-method}
%
\alias{coerce,matrix,ldenseMatrix-method}
\alias{coerce,vector,ldenseMatrix-method}
%
\alias{coerce,matrix,lsparseMatrix-method}
\alias{coerce,vector,lsparseMatrix-method}
%
\description{TODO.}
\examples{
\dontshow{ % for R_DEFAULT_PACKAGES=NULL
library(stats, pos = "package:base", verbose = FALSE)
}
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
