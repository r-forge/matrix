\name{CsparseMatrix-class}
\title{Virtual Class \sQuote{CsparseMatrix} of Compressed Sparse Column Format Matrices}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{CsparseMatrix-class}
%
\alias{ngCMatrix-class}
\alias{lgCMatrix-class}
\alias{igCMatrix-class}
\alias{dgCMatrix-class}
\alias{zgCMatrix-class}
%
\alias{nsCMatrix-class}
\alias{lsCMatrix-class}
\alias{isCMatrix-class}
\alias{dsCMatrix-class}
\alias{zsCMatrix-class}
%
\alias{dpCMatrix-class}
\alias{zpCMatrix-class}
%
\alias{ntCMatrix-class}
\alias{ltCMatrix-class}
\alias{itCMatrix-class}
\alias{dtCMatrix-class}
\alias{ztCMatrix-class}
%
\alias{.validateCsparse}
%
\description{
  \code{CsparseMatrix} is a virtual subclass of
  \code{\linkS4class{sparseMatrix}} representing compressed sparse
  column (CSC) format matrices.  To represent an \code{m}-by-\code{n}
  matrix \code{A}, this format stores the positions of \code{m*n} or
  fewer (often \eqn{O(\max\{m, n\})}{O(max(m, n))}) entries of \code{A},
  ordering these by column then by row.  Authors of subclasses of
  \code{CsparseMatrix} are free to declare how this information is used
  to specify \code{A}.  See \sQuote{Subclasses} for subclasses defined
  in package \pkg{Matrix}.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}}{inherited from virtual superclass
      \code{\linkS4class{Matrix}}.}
    \item{\code{p}}{an integer vector of length \code{Dim[2]+1}
      corresponding to a subset \code{V} of the set of entries of
      \code{A}.  \code{p[j+1]} is the number of elements of \code{V}
      belonging to one of the first \code{j} columns of \code{A},
      implying
      \preformatted{0 = p[1] <= p[j] <= p[j+1] <= p[length(p)] = |V|}
      where \code{|V|} denotes the size of \code{V}.  The number of
      elements of \code{V} belonging to column \code{j} of \code{A} is
      given by \code{p[j+1]-p[j]} or, equivalently, \code{diff(p)[j]}.}
    \item{\code{i}}{an integer vector of length \code{p[length(p)]}
      storing the 0-based row indices of the elements of \code{V}
      ordered by column then by row.  Hence, to be valid, \code{i} must
      satisfy \code{0 <= i < Dim[1]} and be increasing within the
      contiguous segments that partition it by column of \code{A}.
      The segment corresponding to column \code{j} can be accessed as
      \code{i[(p[j]+1):p[j+1]]} if \code{p[j] < p[j+1]} and by
      \code{i[\link{seq.int}(from=p[j], length.out=p[j+1]-p[j])]} more
      generally.  Experimentally, \code{length(i)} is allowed to exceed
      \code{p[length(p)]}; if it does, then the excess elements of
      \code{i} are unreferenced.}
  }
}
\section{Subclasses}{
  The nonvirtual subclasses of \code{CsparseMatrix} defined in package
  \pkg{Matrix} are listed below with links to virtual superclasses
  defining their data type and structure.
  \tabular{rrr}{
    Nonvirtual subclass \tab Data type \tab Structure\cr
    \code{ngCMatrix} \tab \code{\linkS4class{nMatrix}} \tab    \code{\linkS4class{generalMatrix}}\cr
    \code{lgCMatrix} \tab \code{\linkS4class{lMatrix}} \tab    \code{\linkS4class{generalMatrix}}\cr
    \code{igCMatrix} \tab \code{\linkS4class{iMatrix}} \tab    \code{\linkS4class{generalMatrix}}\cr
    \code{dgCMatrix} \tab \code{\linkS4class{dMatrix}} \tab    \code{\linkS4class{generalMatrix}}\cr
    \code{zgCMatrix} \tab \code{\linkS4class{zMatrix}} \tab    \code{\linkS4class{generalMatrix}}\cr
    \code{nsCMatrix} \tab \code{\linkS4class{nMatrix}} \tab  \code{\linkS4class{symmetricMatrix}}\cr
    \code{lsCMatrix} \tab \code{\linkS4class{lMatrix}} \tab  \code{\linkS4class{symmetricMatrix}}\cr
    \code{isCMatrix} \tab \code{\linkS4class{iMatrix}} \tab  \code{\linkS4class{symmetricMatrix}}\cr
    \code{dsCMatrix} \tab \code{\linkS4class{dMatrix}} \tab  \code{\linkS4class{symmetricMatrix}}\cr
    \code{zsCMatrix} \tab \code{\linkS4class{zMatrix}} \tab  \code{\linkS4class{symmetricMatrix}}\cr
    \code{dpCMatrix} \tab \code{\linkS4class{dMatrix}} \tab     \code{\linkS4class{posdefMatrix}}\cr
    \code{zpCMatrix} \tab \code{\linkS4class{zMatrix}} \tab     \code{\linkS4class{posdefMatrix}}\cr
    \code{ntCMatrix} \tab \code{\linkS4class{nMatrix}} \tab \code{\linkS4class{triangularMatrix}}\cr
    \code{ltCMatrix} \tab \code{\linkS4class{lMatrix}} \tab \code{\linkS4class{triangularMatrix}}\cr
    \code{itCMatrix} \tab \code{\linkS4class{iMatrix}} \tab \code{\linkS4class{triangularMatrix}}\cr
    \code{dtCMatrix} \tab \code{\linkS4class{dMatrix}} \tab \code{\linkS4class{triangularMatrix}}\cr
    \code{ztCMatrix} \tab \code{\linkS4class{zMatrix}} \tab \code{\linkS4class{triangularMatrix}}
  }
  The subclasses of \code{iMatrix}, \code{zMatrix}, and
  \code{posdefMatrix} were not defined until \pkg{Matrix} version 1.8-0.

  Classes \code{.gCMatrix} represent general (that is, unstructured)
  \code{A}.  Objects have \code{p} and \code{i} slots specifying a
  subset \code{V} of the set of entries of \code{A}, as described in
  \sQuote{Slots}.  For \code{ngCMatrix}, whose data type is boolean, the
  entries in \code{V} have value 1, and the entries not in \code{V} have
  value 0.  For other \code{.gCMatrix}, the entries in \code{V} have
  values taken from an additional \code{x} slot (a logical, integer,
  double, or complex vector), and the entries not in \code{V} have value
  0.  The \code{i} and \code{x} slots are entirely parallel, both using
  the ordering of \code{V} by column then by row.

  Classes \code{.sCMatrix} represent Hermitian or symmetric \code{A}.
  Compared to \code{.gCMatrix}, objects have an additional \code{uplo}
  slot indicating a restriction on \code{V}.  If \code{uplo="U"}, then
  \code{V} contains only upper triangular entries of \code{A}.  If
  \code{uplo="L"}, then \code{V} contains only lower triangular entries
  of \code{A}.  The entries in \code{V} have value 1 (\code{nsCMatrix})
  or values taken from an \code{x} slot (other \code{.sCMatrix}).  The
  entries opposite the main diagonal, namely \code{A[j, i]} for all
  \code{A[i, j]} in \code{V}, have value \code{op(A[i, j])}, where
  \code{op=Conj} for Hermitian \code{A} and \code{op=identity} for
  symmetric \code{A}.  All other entries of \code{A} have value 0.
  Objects of class \code{zsCMatrix} have an additional \code{trans}
  slot that distinguishes Hermitian and symmetric \code{A}, which use
  \code{trans="C"} and \code{"T"}, respectively.  If \code{trans="C"}
  and \code{V} contains diagonal entries of \code{A}, then the imaginary
  parts of the corresponding elements of \code{x} are not referenced.
  (The diagonal entries of Hermitian matrices are real by definition.)

  Classes \code{dpCMatrix} and \code{zpCMatrix} extend \code{dsCMatrix}
  and \code{zsCMatrix} and represent positive semidefinite \code{A}.
  Positive semidefiniteness of \code{A} induces a constraint on the
  \code{x} slot.  However, this constraint is not strictly enforced by
  the validity methods, which test only that the diagonal entries of
  \code{A} are nonnegative (a necessary, not sufficient condition).
  For \code{zpCMatrix}, only \code{trans="C"} is valid.

  Classes \code{.tCMatrix} represent triangular \code{A}.  Compared to
  \code{.gCMatrix}, objects have additional \code{uplo} and \code{diag}
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
  \dQuote{Sibling} classes \code{\linkS4class{RsparseMatrix}} and
  \code{\linkS4class{TsparseMatrix}} for compressed sparse row and
  triplet format matrices.  Generic functions
  \code{\linkS4methods{cbind2}}, \code{\linkS4methods{colSums}},
  \code{\linkS4methods{[}} (when called as \code{x[, j]}), etc., whose
  methods for \code{CsparseMatrix} are particularly efficient.
}
\examples{
getClass("CsparseMatrix")

## The common validity check function (based on C code):
getValidity(getClass("CsparseMatrix"))

%% moved from ./dgCMatrix-class.Rd :
(m <- Matrix(c(0,0,2:0), 3,5))
str(m)
m[,1]
\dontshow{## regression test: this must give a validity-check error:
stopifnot(inherits(try(new("dgCMatrix", i = 0:1, p = 0:2,
                           x = c(2,3), Dim = 3:4)),
          "try-error"))
}

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
