\name{TsparseMatrix-class}
\title{Virtual Class \sQuote{TsparseMatrix} of Triplet Format Matrices}
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
\alias{dpTMatrix-class}
\alias{zpTMatrix-class}
%
\alias{ntTMatrix-class}
\alias{ltTMatrix-class}
\alias{itTMatrix-class}
\alias{dtTMatrix-class}
\alias{ztTMatrix-class}
%
\description{
  \code{TsparseMatrix} is a virtual subclass of
  \code{\linkS4class{sparseMatrix}} representing triplet format
  matrices.  To represent an \code{m}-by-\code{n} matrix \code{A}, this
  format stores the positions of \code{m*n} or fewer
  (often \eqn{O(\max\{m, n\})}{O(max(m, n))}) entries of \code{A}.
  Authors of subclasses of \code{TsparseMatrix} are free to declare how
  this information is used to specify \code{A}.  See \sQuote{Subclasses}
  for subclasses defined in package \pkg{Matrix}.

  The triplet format is convenient for reading and writing sparse
  matrices but often not amenable to direct computation.  Accordingly,
  methods for \code{TsparseMatrix} defined in package \pkg{Matrix} tend
  to work by coercion to and from a more tractable \emph{compressed}
  sparse format.  Users can avoid repeated coercions by methods by
  coercing sparse matrix operands to either
  \code{\linkS4class{CsparseMatrix}} or
  \code{\linkS4class{RsparseMatrix}} \emph{before} computing with them.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}}{inherited from virtual superclass
      \code{\linkS4class{Matrix}}.}
    \item{\code{i}, \code{j}}{integer vectors of equal length
      corresponding to a subset \code{V} of the set of entries of
      \code{A}.  Respectively, \code{i} and \code{j} store the 0-based
      row and column indices of the elements of \code{V}.  Hence, to be
      valid, \code{i} and \code{j} must satisfy \code{0 <= i < Dim[1]}
      and \code{0 <= j < Dim[2]}.  The common length of \code{i} and
      \code{j} can exceed the size of \code{V}, allowing for duplication
      of the form \code{i[k1] == i[k2]}, \code{j[k1] == j[k2]},
      \code{k1 != k2}.}
  }
}
\section{Subclasses}{
  The nonvirtual subclasses of \code{TsparseMatrix} defined in package
  \pkg{Matrix} are listed below with links to virtual superclasses
  defining their data type and structure.
  \tabular{rrr}{
    Nonvirtual subclass \tab Data type \tab Structure\cr
    \code{ngTMatrix} \tab \code{\linkS4class{nMatrix}} \tab    \code{\linkS4class{generalMatrix}}\cr
    \code{lgTMatrix} \tab \code{\linkS4class{lMatrix}} \tab    \code{\linkS4class{generalMatrix}}\cr
    \code{igTMatrix} \tab \code{\linkS4class{iMatrix}} \tab    \code{\linkS4class{generalMatrix}}\cr
    \code{dgTMatrix} \tab \code{\linkS4class{dMatrix}} \tab    \code{\linkS4class{generalMatrix}}\cr
    \code{zgTMatrix} \tab \code{\linkS4class{zMatrix}} \tab    \code{\linkS4class{generalMatrix}}\cr
    \code{nsTMatrix} \tab \code{\linkS4class{nMatrix}} \tab  \code{\linkS4class{symmetricMatrix}}\cr
    \code{lsTMatrix} \tab \code{\linkS4class{lMatrix}} \tab  \code{\linkS4class{symmetricMatrix}}\cr
    \code{isTMatrix} \tab \code{\linkS4class{iMatrix}} \tab  \code{\linkS4class{symmetricMatrix}}\cr
    \code{dsTMatrix} \tab \code{\linkS4class{dMatrix}} \tab  \code{\linkS4class{symmetricMatrix}}\cr
    \code{zsTMatrix} \tab \code{\linkS4class{zMatrix}} \tab  \code{\linkS4class{symmetricMatrix}}\cr
    \code{dpTMatrix} \tab \code{\linkS4class{dMatrix}} \tab     \code{\linkS4class{posdefMatrix}}\cr
    \code{zpTMatrix} \tab \code{\linkS4class{zMatrix}} \tab     \code{\linkS4class{posdefMatrix}}\cr
    \code{ntTMatrix} \tab \code{\linkS4class{nMatrix}} \tab \code{\linkS4class{triangularMatrix}}\cr
    \code{ltTMatrix} \tab \code{\linkS4class{lMatrix}} \tab \code{\linkS4class{triangularMatrix}}\cr
    \code{itTMatrix} \tab \code{\linkS4class{iMatrix}} \tab \code{\linkS4class{triangularMatrix}}\cr
    \code{dtTMatrix} \tab \code{\linkS4class{dMatrix}} \tab \code{\linkS4class{triangularMatrix}}\cr
    \code{ztTMatrix} \tab \code{\linkS4class{zMatrix}} \tab \code{\linkS4class{triangularMatrix}}
  }
  The subclasses of \code{iMatrix}, \code{zMatrix}, and
  \code{posdefMatrix} were not defined until \pkg{Matrix} version 1.8-0.

  Classes \code{.gTMatrix} represent general (that is, unstructured)
  \code{A}.  Objects have \code{i} and \code{j} slots specifying a
  subset \code{V} of the set of entries of \code{A}, as described in
  \sQuote{Slots}.  For \code{ngTMatrix}, whose data type is boolean, the
  entries in \code{V} have value 1, and the entries not in \code{V} have
  value 0.  For other \code{.gTMatrix}, the entries in \code{V} have
  values taken from an additional \code{x} slot (a logical, integer,
  double, or complex vector), and the entries not in \code{V} have value
  0.  The \code{i}, \code{j}, and \code{x} slots have a common length,
  and the value of \code{A[i[k], j[k]]} is the result of aggregating
  \code{x[i == i[k] & j == j[k]]} by summation (reduction by \code{|} if
  \code{x} is logical, by \code{+} otherwise).

  Classes \code{.sCMatrix} represent Hermitian or symmetric \code{A}.
  Compared to \code{.gTMatrix}, objects have an additional \code{uplo}
  slot indicating a restriction on \code{V}.  If \code{uplo="U"}, then
  \code{V} contains only upper triangular entries of \code{A}.  If
  \code{uplo="L"}, then \code{V} contains only lower triangular entries
  of \code{A}.  The entries in \code{V} have value 1 (\code{nsCMatrix})
  or values taken from an \code{x} slot (other \code{.sCMatrix}).  The
  entries opposite the main diagonal, namely \code{A[j, i]} for all
  \code{A[i, j]} in \code{V}, have value \code{op(A[i, j])}, where
  \code{op=Conj} for Hermitian \code{A} and \code{op=identity} for
  symmetric \code{A}.  All other entries of \code{A} have value 0.
  Objects of class \code{zsTMatrix} have an additional \code{trans}
  slot that distinguishes Hermitian and symmetric \code{A}, which use
  \code{trans="C"} and \code{"T"}, respectively.  If \code{trans="C"}
  and \code{V} contains diagonal entries of \code{A}, then the imaginary
  parts of the corresponding elements of \code{x} are not referenced.
  (The diagonal entries of Hermitian matrices are real by definition.)

  Classes \code{dpTMatrix} and \code{zpTMatrix} extend \code{dsTMatrix}
  and \code{zsTMatrix} and represent positive semidefinite \code{A}.
  Positive semidefiniteness of \code{A} induces a constraint on the
  \code{x} slot.  However, this constraint is not strictly enforced by
  the validity methods, which test only that the diagonal entries of
  \code{A} are nonnegative (a necessary, not sufficient condition).
  For \code{zpTMatrix}, only \code{trans="C"} is valid.

  Classes \code{.tTMatrix} represent triangular \code{A}.  Compared to
  \code{.gTMatrix}, objects have additional \code{uplo} and \code{diag}
  slots indicating restrictions on \code{V}.  If \code{uplo="U"}, then
  \code{A} is upper triangular and \code{V} contains only upper
  triangular entries.  If \code{uplo="L"}, then \code{A} is lower
  triangular and \code{V} contains only lower triangular entries.  If
  \code{diag="U"} (instead of \code{"N"}), then \code{A} is unit
  triangular and \code{V} contains only \emph{strictly} upper or lower
  triangular entries.  The entries in \code{V} have value 1
  (\code{ntCMatrix}) or values taken from an \code{x} slot (other
  \code{.tCMatrix}).  With the exception of diagonal entries when
  \code{diag="U"}, which have value 1, entries not in \code{V} have
  value 0.
}
\seealso{
  \dQuote{Sibling} classes \code{\linkS4class{CsparseMatrix}} and
  \code{\linkS4class{RsparseMatrix}} for compressed sparse column and
  compressed sparse row format matrices.  Functions
  \code{\link{aggregateT}} and \code{\link{asUniqueT}} for
  \dQuote{canonicalizing} objects inheriting from \code{TsparseMatrix}.
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
