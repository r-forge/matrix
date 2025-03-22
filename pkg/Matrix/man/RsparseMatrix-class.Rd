\name{RsparseMatrix-class}
\title{Virtual Class \sQuote{RsparseMatrix} of Compressed Sparse Row Format Matrices}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{RsparseMatrix-class}
%
\alias{ngRMatrix-class}
\alias{lgRMatrix-class}
\alias{igRMatrix-class}
\alias{dgRMatrix-class}
\alias{zgRMatrix-class}
%
\alias{nsRMatrix-class}
\alias{lsRMatrix-class}
\alias{isRMatrix-class}
\alias{dsRMatrix-class}
\alias{zsRMatrix-class}
%
\alias{dpRMatrix-class}
\alias{zpRMatrix-class}
%
\alias{ntRMatrix-class}
\alias{ltRMatrix-class}
\alias{itRMatrix-class}
\alias{dtRMatrix-class}
\alias{ztRMatrix-class}
%
\description{
  \code{RsparseMatrix} is a virtual subclass of
  \code{\linkS4class{sparseMatrix}} representing compressed sparse row
  (CSR) format matrices.  To represent an \code{m}-by-\code{n} matrix
  \code{A}, this format stores the positions of \code{m*n} or fewer
  (often \eqn{O(\max\{m, n\})}{O(max(m, n))}) entries of \code{A},
  ordering these by row then by column.  Authors of subclasses of
  \code{RsparseMatrix} are free to declare how this information is used
  to specify \code{A}.  See \sQuote{Subclasses} for subclasses defined
  in package \pkg{Matrix}.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}}{inherited from virtual superclass
      \code{\linkS4class{Matrix}}.}
    \item{\code{p}}{an integer vector of length \code{Dim[1]+1}
      corresponding to a subset \code{V} of the set of entries of
      \code{A}.  \code{p[i+1]} is the number of elements of \code{V}
      belonging to one of the first \code{i} rows of \code{A}, implying
      \preformatted{0 = p[1] <= p[i] <= p[i+1] <= p[length(p)] = |V|}
      where \code{|V|} denotes the size of \code{V}.  The number of
      elements of \code{V} belonging to row \code{i} of \code{A} is
      given by \code{p[i+1]-p[i]} or, equivalently, \code{diff(p)[i]}.}
    \item{\code{j}}{an integer vector of length \code{p[length(p)]}
      storing the 0-based column indices of the elements of \code{V}
      ordered by row then by column.  Hence, to be valid, \code{j} must
      satisfy \code{0 <= j < Dim[2]} and be increasing within the
      contiguous segments that partition it by row of \code{A}.  The
      segment corresponding to row \code{i} can be accessed as
      \code{j[(p[i]+1):p[i+1]]} if \code{p[i] < p[i+1]} and by
      \code{j[\link{seq.int}(from=p[i], length.out=p[i+1]-p[i])]} more
      generally.  Experimentally, \code{length(j)} is allowed to exceed
      \code{p[length(p)]}; if it does, then the excess elements of
      \code{j} are unreferenced.}
  }
}
\section{Subclasses}{
  The nonvirtual subclasses of \code{RsparseMatrix} defined in package
  \pkg{Matrix} are listed below with links to virtual superclasses
  defining their data type and structure.
  \tabular{rrr}{
    Nonvirtual subclass \tab Data type \tab Structure\cr
    \code{ngRMatrix} \tab \code{\linkS4class{nMatrix}} \tab    \code{\linkS4class{generalMatrix}}\cr
    \code{lgRMatrix} \tab \code{\linkS4class{lMatrix}} \tab    \code{\linkS4class{generalMatrix}}\cr
    \code{igRMatrix} \tab \code{\linkS4class{iMatrix}} \tab    \code{\linkS4class{generalMatrix}}\cr
    \code{dgRMatrix} \tab \code{\linkS4class{dMatrix}} \tab    \code{\linkS4class{generalMatrix}}\cr
    \code{zgRMatrix} \tab \code{\linkS4class{zMatrix}} \tab    \code{\linkS4class{generalMatrix}}\cr
    \code{nsRMatrix} \tab \code{\linkS4class{nMatrix}} \tab  \code{\linkS4class{symmetricMatrix}}\cr
    \code{lsRMatrix} \tab \code{\linkS4class{lMatrix}} \tab  \code{\linkS4class{symmetricMatrix}}\cr
    \code{isRMatrix} \tab \code{\linkS4class{iMatrix}} \tab  \code{\linkS4class{symmetricMatrix}}\cr
    \code{dsRMatrix} \tab \code{\linkS4class{dMatrix}} \tab  \code{\linkS4class{symmetricMatrix}}\cr
    \code{zsRMatrix} \tab \code{\linkS4class{zMatrix}} \tab  \code{\linkS4class{symmetricMatrix}}\cr
    \code{dpRMatrix} \tab \code{\linkS4class{dMatrix}} \tab     \code{\linkS4class{posdefMatrix}}\cr
    \code{zpRMatrix} \tab \code{\linkS4class{zMatrix}} \tab     \code{\linkS4class{posdefMatrix}}\cr
    \code{ntRMatrix} \tab \code{\linkS4class{nMatrix}} \tab \code{\linkS4class{triangularMatrix}}\cr
    \code{ltRMatrix} \tab \code{\linkS4class{lMatrix}} \tab \code{\linkS4class{triangularMatrix}}\cr
    \code{itRMatrix} \tab \code{\linkS4class{iMatrix}} \tab \code{\linkS4class{triangularMatrix}}\cr
    \code{dtRMatrix} \tab \code{\linkS4class{dMatrix}} \tab \code{\linkS4class{triangularMatrix}}\cr
    \code{ztRMatrix} \tab \code{\linkS4class{zMatrix}} \tab \code{\linkS4class{triangularMatrix}}
  }
  The subclasses of \code{iMatrix}, \code{zMatrix}, and
  \code{posdefMatrix} were not defined until \pkg{Matrix} version 1.8-0.

  Classes \code{.gRMatrix} represent general (unstructured, possibly
  nonsquare) \code{A}.  Objects have \code{p} and \code{j} slots
  specifying a subset \code{V} of the set of entries of \code{A}, as
  described in \sQuote{Slots}.  For \code{ngRMatrix}, whose data type is
  boolean, the entries in \code{V} have value 1, and the entries not in
  \code{V} have value 0.  For other \code{.gRMatrix}, the entries in
  \code{V} have values taken from an additional \code{x} slot (a
  logical, integer, double, or complex vector), and the entries not in
  \code{V} have value 0.  The \code{j} and \code{x} slots are entirely
  parallel, both using the ordering of \code{V} by row then by column.

  Classes \code{.sRMatrix} represent Hermitian or symmetric \code{A}.
  Compared to \code{.gRMatrix}, objects have an additional \code{uplo}
  slot indicating a restriction on \code{V}.  If \code{uplo="U"}, then
  \code{V} contains only upper triangular entries of \code{A}.  If
  \code{uplo="L"}, then \code{V} contains only lower triangular entries
  of \code{A}.  The entries in \code{V} have value 1 (\code{nsRMatrix})
  or values taken from an \code{x} slot (other \code{.sRMatrix}).  The
  entries opposite the main diagonal, namely \code{A[j, i]} for all
  \code{A[i, j]} in \code{V}, have value \code{op(A[i, j])}, where
  \code{op=Conj} for Hermitian \code{A} and \code{op=identity} for
  symmetric \code{A}.  All other entries of \code{A} have value 0.
  Objects of class \code{zsRMatrix} have an additional \code{trans}
  slot that distinguishes Hermitian and symmetric \code{A}, which use
  \code{trans="C"} and \code{"T"}, respectively.  If \code{trans="C"}
  and \code{V} contains diagonal entries of \code{A}, then the imaginary
  parts of the corresponding elements of \code{x} are not referenced.
  (The diagonal entries of Hermitian matrices are real by definition.)

  Classes \code{dpRMatrix} and \code{zpRMatrix} extend \code{dsRMatrix}
  and \code{zsRMatrix} and represent positive semidefinite \code{A}.
  Positive semidefiniteness of \code{A} induces a constraint on the
  \code{x} slot.  However, this constraint is not strictly enforced by
  the validity methods, which test only that the diagonal entries of
  \code{A} are nonnegative (a necessary, not sufficient condition).
  For \code{zpRMatrix}, only \code{trans="C"} is valid.

  Classes \code{.tRMatrix} represent triangular \code{A}.  Compared to
  \code{.gRMatrix}, objects have additional \code{uplo} and \code{diag}
  slots indicating restrictions on \code{V}.  If \code{uplo="U"}, then
  \code{A} is upper triangular and \code{V} contains only upper
  triangular entries.  If \code{uplo="L"}, then \code{A} is lower
  triangular and \code{V} contains only lower triangular entries.  If
  \code{diag="U"} (instead of \code{"N"}), then \code{A} is unit
  triangular and \code{V} contains only \emph{strictly} upper or lower
  triangular entries.  With the exception of diagonal entries when
  \code{diag="U"}, which have value 1, entries not in \code{V} have
  value 0.
}
\seealso{
  \dQuote{Sister} classes \code{\linkS4class{CsparseMatrix}} and
  \code{\linkS4class{TsparseMatrix}} for compressed sparse column and
  triplet format matrices.  Generic functions
  \code{\linkS4methods{rbind2}}, \code{\linkS4methods{rowSums}},
  \code{\linkS4methods{[}} (when called as \code{x[i, ]}), etc., whose
  methods for \code{RsparseMatrix} are particularly efficient.
}
\examples{
showClass("RsparseMatrix")

%% moved from ./dsRMatrix-class.Rd :
(m0 <- new("dsRMatrix"))
m2 <- new("dsRMatrix", Dim = c(2L,2L),
          x = c(3,1), j = c(1L,1L), p = 0:2)
m2
stopifnot(colSums(as(m2, "TsparseMatrix")) == 3:4)
str(m2)
(ds2 <- forceSymmetric(diag(2))) # dsy*
dR <- as(ds2, "RsparseMatrix")
dR # dsRMatrix

%% moved from ./dtRMatrix-class.Rd :
(m0 <- new("dtRMatrix"))
(m2 <- new("dtRMatrix", Dim = c(2L,2L),
                        x = c(5, 1:2), p = c(0L,2:3), j= c(0:1,1L)))
str(m2)
(m3 <- as(Diagonal(2), "RsparseMatrix"))# --> dtRMatrix
}
