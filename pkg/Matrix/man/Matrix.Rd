\name{Matrix}
\title{Construct a Classed Matrix}
%
\keyword{array}
\keyword{utilities}
%
\alias{Matrix}
\description{
  Construct a Matrix of a class inheriting from \code{Matrix},
  i.e., a \dQuote{\pkg{Matrix}-package \code{Matrix} object}, using a
  syntax known from base \R's \code{\link{matrix}()}.
}
\usage{
Matrix(data=NA, nrow=1, ncol=1, byrow=FALSE, dimnames=NULL,
       sparse = NULL, doDiag = TRUE, forceCheck = FALSE,
       int2dbl = !as.logical(getOption("Matrix.use_iMatrix",
                                       .MatrixEnv[["use_iMatrix"]])))
% future default: int2dbl = FALSE	%
}
\arguments{
  \item{data}{an optional \dQuote{numeric-alike} data vector,
    \code{\link{matrix}} (or also \code{"\linkS4class{Matrix}"}).}
  \item{nrow, ncol}{when \code{data} is not a matrix, the desired number of
    rows and columns, respectively.}
  \item{byrow}{logical.  If \code{FALSE} (the default) the matrix is
    filled by columns, otherwise the matrix is filled by rows.}
  \item{dimnames}{a \code{\link{dimnames}} attribute for the matrix: a
    \code{list} of two character components.  They are set if not
    \code{\link{NULL}} (as per default).}
  \item{sparse}{logical or \code{NULL}, specifying if the result should
    be sparse or not.  By default, it is made sparse when more than half
    of the entries are 0.}
  \item{doDiag}{logical indicating if a \code{\linkS4class{diagonalMatrix}}
    object should be returned when the resulting matrix is diagonal
    (\emph{mathematically}).  As class \code{\linkS4class{diagonalMatrix}}
    \code{\link{extends}} \code{\linkS4class{sparseMatrix}}, this is a
    natural default for all values of \code{sparse}.

    Otherwise, if \code{doDiag} is false, a dense or sparse (depending on
    \code{sparse}) \emph{symmetric} matrix will be returned.}
% MJ: not quite ... if 'dimnames' is asymmetric then the result
%     is a triangularMatrix, not a symmetricMatrix ...

  \item{forceCheck}{logical indicating if the checks for structure
    should even happen when \code{data} is already a \code{"Matrix"}
    object.}
  \item{int2dbl}{optional \code{\link{logical}} argument (new in
    \pkg{Matrix} 1.8.0), determining if \code{\link{integer}} \code{data}
    should be treated as \code{\link{double}}, back-compatibly to the pre-1.8.0 
    behaviour, but not desirable any longer with \code{"\linkS4class{iMatrix}"}
    (integer valued) Matrix classes available now.
    While the default currently depends on the \code{\link{options}("Matrix.use_Matrix_i")}
    setting, the default will become \code{FALSE} in the future.}
}
\value{
  Returns matrix of a class that inherits from \code{"\linkS4class{Matrix}"}.
  Only if \code{data} is not a \code{\link{matrix}} and does not already inherit
  from class \code{\linkS4class{Matrix}} are the arguments
  \code{nrow}, \code{ncol} and \code{byrow} made use of.
}
\details{
  If either of \code{nrow} or \code{ncol} is not given, an attempt is
  made to infer it from the length of \code{data} and the other
  parameter (analogously to how \code{\link{matrix}()} works).
  Further, \code{Matrix()} makes efforts to keep \code{\link{logical}}
  matrices logical, i.e., inheriting from class \code{\linkS4class{lMatrix}},
  and analogously for \code{integer} or \code{complex} \code{data}.
  \code{\link{raw}} \code{data} is currently coerced to
  \code{\link{integer}} and treated as such.
  \code{Matrix()} also determines specially structured matrices such as symmetric,
  triangular or diagonal ones.  Note that a \emph{symmetric} matrix also
  needs symmetric \code{\link{dimnames}}, e.g., by specifying
  \code{dimnames = list(NULL,NULL)}, see the examples.

  Most of the time, the function works via a traditional (\emph{full})
  \code{\link{matrix}}.  However, \code{Matrix(0, nrow,ncol)} directly
  constructs an \dQuote{empty} \linkS4class{sparseMatrix}, as does
  \code{Matrix(FALSE, *)}.

  Although it is sometime possible to mix unclassed matrices (created
  with \code{matrix}) with ones of class \code{"Matrix"}, it is much
  safer to always use carefully constructed ones of class
  \code{"\linkS4class{Matrix}"}.
}
\seealso{
  The classes \code{\linkS4class{Matrix}},
  \code{\linkS4class{symmetricMatrix}},
  \code{\linkS4class{triangularMatrix}}, and
  \code{\linkS4class{diagonalMatrix}}; further,
  \code{\link{matrix}}.

  Special matrices can be constructed, e.g., via
  \code{\link{sparseMatrix}} (sparse), \code{\link{bdiag}}
  (block-diagonal), \code{\link{bandSparse}} (banded sparse), or
  \code{\link{Diagonal}}.
}
\examples{
\dontshow{ % for R_DEFAULT_PACKAGES=NULL %
library(utils, pos = "package:base", verbose = FALSE)
}
Matrix(0, 3, 2)             # 3 by 2 matrix of zeros -> sparse
Matrix(0, 3, 2, sparse=FALSE)# -> 'dense'

## 4 cases - 3 different results :
Matrix(0, 2, 2)              # diagonal !
Matrix(0, 2, 2, sparse=FALSE)# (ditto)
Matrix(0, 2, 2,               doDiag=FALSE)# -> sparse symm. "dsCMatrix"
Matrix(0, 2, 2, sparse=FALSE, doDiag=FALSE)# -> dense  symm. "dsyMatrix"

Matrix(1:6, 3, 2)           # a 3 by 2 matrix (+ integer warning)
Matrix(1:6 + 1, nrow=3)

## logical ones:
Matrix(diag(4) >  0) # -> "ldiMatrix" with diag = "U"
Matrix(diag(4) >  0, sparse=TRUE) #  (ditto)
Matrix(diag(4) >= 0) # -> "lsyMatrix" (of all 'TRUE')
## triangular
l3 <- upper.tri(matrix(,3,3))
(M <- Matrix(l3))   # -> "ltCMatrix"
Matrix(! l3)        # -> "ltrMatrix"
as(l3, "CsparseMatrix")# "lgCMatrix"

Matrix(1:9, nrow=3,
       dimnames = list(c("a", "b", "c"), c("A", "B", "C")))
(I3 <- Matrix(diag(3)))# identity, i.e., unit "diagonalMatrix"
str(I3) # note  'diag = "U"' and the empty 'x' slot

(A <- cbind(a=c(2,1), b=1:2))# symmetric *apart* from dimnames
Matrix(A)                    # hence 'dgeMatrix'
(As <- Matrix(A, dimnames = list(NULL,NULL)))# -> symmetric
forceSymmetric(A) # also symmetric, w/ symm. dimnames
stopifnot(is(As, "symmetricMatrix"),
          is(Matrix(0, 3,3), "sparseMatrix"),
          is(Matrix(FALSE, 1,1), "sparseMatrix"))
}
