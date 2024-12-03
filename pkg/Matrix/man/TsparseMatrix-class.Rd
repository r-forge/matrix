\name{TsparseMatrix-class}
\title{Class "TsparseMatrix" of Sparse Matrices in Triplet Form}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{TsparseMatrix-class}
%
\alias{ngTMatrix-class}
\alias{lgTMatrix-class}
\alias{igTMatrix-class}
\alias{dgTMatrix-class}
\alias{zgTMatrix-class}
%
\alias{nsTMatrix-class}
\alias{lsTMatrix-class}
\alias{isTMatrix-class}
\alias{dsTMatrix-class}
\alias{zsTMatrix-class}
%
\alias{ntTMatrix-class}
\alias{ltTMatrix-class}
\alias{itTMatrix-class}
\alias{dtTMatrix-class}
\alias{ztTMatrix-class}
%
\description{The \code{"TsparseMatrix"} class is the virtual class of
  all sparse matrices coded in triplet form.  Since it is a virtual class,
  no objects may be created from it.  See
  \code{showClass("TsparseMatrix")} for its subclasses.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}:}{from the \code{"\linkS4class{Matrix}"} class,}
    \item{\code{i}:}{Object of class \code{"integer"} - the row indices
      of non-zero entries \emph{in 0-base}, i.e., must be in
      \code{0:(nrow(.)-1)}.}
    \item{\code{j}:}{Object of class \code{"integer"} - the column
      indices of non-zero entries.  Must be the same length as slot
      \code{i} and \emph{0-based} as well, i.e., in
      \code{0:(ncol(.)-1)}.   For numeric Tsparse matrices, \code{(i,j)}
      pairs can occur more than once, see \code{\linkS4class{dgTMatrix}}.
    }
  }
}
\section{Extends}{
Class \code{"sparseMatrix"}, directly.
Class \code{"Matrix"}, by class \code{"sparseMatrix"}.
}
\section{Methods}{
  Extraction (\code{"["}) methods, see
  \code{\link{[-methods}}.%-> ./Xtrct-methods.Rd
}
\note{
  Most operations with sparse matrices are performed using the
  compressed, column-oriented or \code{\linkS4class{CsparseMatrix}}
  representation.  The triplet representation is convenient for
  creating a sparse matrix or for reading and writing such
  matrices.  Once it is created, however, the matrix is generally
  coerced to a \code{\linkS4class{CsparseMatrix}} for further
  operations.

  Note that all \code{new(.)}, \code{\link{spMatrix}} and
  \code{\link{sparseMatrix}(*, repr="T")} constructors
  for \code{"TsparseMatrix"} classes implicitly add (i.e., \dQuote{sum up})
  \eqn{x_k}'s that belong to identical \eqn{(i_k, j_k)} pairs, see, the
  example below, or also \code{"\linkS4class{dgTMatrix}"}.

  For convenience, methods for some operations such as \code{\%*\%}
  and \code{crossprod} are defined for
  \code{\linkS4class{TsparseMatrix}} objects.  These methods simply
  coerce the \code{\linkS4class{TsparseMatrix}} object to a
  \code{\linkS4class{CsparseMatrix}} object then perform the
  operation.
}
% \author{Martin Maechler}
\seealso{
  its superclass, \code{\linkS4class{sparseMatrix}}, and the
  \code{\linkS4class{dgTMatrix}} class, for the links to other classes.
}
\examples{
showClass("TsparseMatrix")
## or just the subclasses' names
names(getClass("TsparseMatrix")@subclasses)

T3 <- spMatrix(3,4, i=c(1,3:1), j=c(2,4:2), x=1:4)
T3 # only 3 non-zero entries, 5 = 1+4 !
\dontshow{stopifnot(nnzero(T3) == 3)}

%% moved from ./dgTMatrix-class.Rd :
m <- Matrix(0+1:28, nrow = 4)
m[-3,c(2,4:5,7)] <- m[ 3, 1:4] <- m[1:3, 6] <- 0
(mT <- as(m, "TsparseMatrix"))
str(mT)
mT[1,]
mT[4, drop = FALSE]
stopifnot(identical(mT[lower.tri(mT)],
                    m [lower.tri(m) ]))
mT[lower.tri(mT,diag=TRUE)] <- 0
mT

## Triplet representation with repeated (i,j) entries
## *adds* the corresponding x's:
T2 <- new("dgTMatrix",
          i = as.integer(c(1,1,0,3,3)),
          j = as.integer(c(2,2,4,0,0)), x=10*1:5, Dim=4:5)
str(T2) # contains (i,j,x) slots exactly as above, but
T2 ## has only three non-zero entries, as for repeated (i,j)'s,
   ## the corresponding x's are "implicitly" added
stopifnot(nnzero(T2) == 3)

%% moved from ./dsCMatrix-class.Rd :
mm <- Matrix(toeplitz(c(10, 0, 1, 0, 3)), sparse = TRUE)
mm # automatically dsCMatrix
str(mm)
mT <- as(as(mm, "generalMatrix"), "TsparseMatrix")

## Either
(symM <- as(mT, "symmetricMatrix")) # dsT
(symC <- as(symM, "CsparseMatrix")) # dsC
## or
sT <- Matrix(mT, sparse=TRUE, forceCheck=TRUE) # dsT

sym2 <- as(symC, "TsparseMatrix")
## --> the same as 'symM', a "dsTMatrix"
\dontshow{
stopifnot(identical(sT, symM), identical(sym2, symM),
          class(sym2) == "dsTMatrix",
	  identical(sym2[1,], sT[1,]),
	  identical(sym2[,2], sT[,2]))
}

%% moved from ./dtCMatrix-class.Rd :
showClass("dtCMatrix")
showClass("dtTMatrix")
t1 <- new("dtTMatrix", x= c(3,7), i= 0:1, j=3:2, Dim= as.integer(c(4,4)))
t1
## from  0-diagonal to unit-diagonal {low-level step}:
tu <- t1 ; tu@diag <- "U"
tu
(cu <- as(tu, "CsparseMatrix"))
str(cu)# only two entries in @i and @x
stopifnot(cu@i == 1:0,
          all(2 * symmpart(cu) == Diagonal(4) + forceSymmetric(cu)))

t1[1,2:3] <- -1:-2
diag(t1) <- 10*c(1:2,3:2)
t1 # still triangular
(it1 <- solve(t1))
t1. <- solve(it1)
all(abs(t1 - t1.) < 10 * .Machine$double.eps)

## 2nd example
U5 <- new("dtCMatrix", i= c(1L, 0:3), p=c(0L,0L,0:2, 5L), Dim = c(5L, 5L),
          x = rep(1, 5), diag = "U")
U5
(iu <- solve(U5)) # contains one '0'
validObject(iu2 <- solve(U5, Diagonal(5)))# failed in earlier versions

I5 <- iu  \%*\% U5 # should equal the identity matrix
i5 <- iu2 \%*\% U5
m53 <- matrix(1:15, 5,3, dimnames=list(NULL,letters[1:3]))
asDiag <- function(M) as(drop0(M), "diagonalMatrix")
stopifnot(
   all.equal(Diagonal(5), asDiag(I5), tolerance=1e-14) ,
   all.equal(Diagonal(5), asDiag(i5), tolerance=1e-14) ,
   identical(list(NULL, dimnames(m53)[[2]]), dimnames(solve(U5, m53)))
)
\dontshow{
i5. <- I5; colnames(i5.) <- LETTERS[11:15]
M53 <- as(m53, "denseMatrix")
stopifnot(
   identical((dns <- dimnames(solve(i5., M53))),
             dimnames(solve(as.matrix(i5.), as.matrix(M53)))) ,
   identical(dns, dimnames(solve(i5., as.matrix(M53))))
)
}
}
