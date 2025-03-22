\name{indexMatrix-class}
\title{Virtual Class \sQuote{indexMatrix} of Generalized Index Matrices}
%
\docType{class}
\keyword{array}
\keyword{classes}
%
\alias{indexMatrix-class}
%
\alias{ninMatrix-class}
\alias{linMatrix-class}
\alias{iinMatrix-class}
\alias{dinMatrix-class}
\alias{zinMatrix-class}
%
\alias{indMatrix-class}
\alias{pMatrix-class}
%
\description{
  \code{indexMatrix} is a virtual subclass of
  \code{\linkS4class{Matrix}} representing generalized row or column
  index matrices.  A generalized row index matrix is a matrix with at
  most one nonzero entry in each row.  A generalized column index matrix
  is a matrix with at most one nonzero entry in each column.

  Multiplying a matrix on the left by a generalized row index matrix is
  equivalent to sampling rows with replacement and scaling the sampled
  rows.  Multiplying a matrix on the right by a generalized column index
  matrix is equivalent to sampling columns with replacement and scaling
  the sampled columns.

  Standard index matrices are generalized index matrices whose rows or
  columns are all standard unit vectors.  Permutation matrices are
  standard index matrices whose rows \emph{and} columns are all standard
  unit vectors.  These special cases are represented by nonvirtual
  subclasses \code{indMatrix} and \code{pMatrix}, documented below.

  Diagonal matrices are generalized index matrices which are square and
  whose nonzero entries are situated on the main diagonal.  This special
  case is represented by virtual subclass
  \code{\linkS4class{diagonalMatrix}}.
}
\section{Slots}{
  \describe{
    \item{\code{Dim}, \code{Dimnames}}{inherited from virtual superclass
      \code{\linkS4class{Matrix}}.}
  }
  All direct, nonvirtual subclasses of \code{indexMatrix} defined in
  package \pkg{Matrix} have two additional slots:
  \describe{
    \item{\code{margin}}{an integer, either 1 or 2, indicating that the
      represented matrix is a generalized row (1) or column (2) index
      matrix.}
    \item{\code{perm}}{an integer vector of length \code{Dim[margin]}
      with elements sampled with replacement from
      \code{seq_len(Dim[-margin])}.  The represented matrix is zero
      except for entries at \code{[i, perm[i]]} (\code{margin=1}) or
      \code{[perm[j], j]} (\code{margin=2}), which may be nonzero.}
  }
}
\section{Subclasses}{
  Currently, the only direct, nonvirtual subclass of \code{indexMatrix}
  is \code{indMatrix},  which represents standard index matrices.
  Objects have \code{margin} and \code{perm} slots, as described in
  \sQuote{Slots}.  The represented matrix is equal to 1 in positions
  \code{[i, perm[i]]} (\code{margin=1}) or \code{[perm[j], j]}
  (\code{margin=2}) and otherwise equal to 0.  Thus, \code{indMatrix}
  extends not only \code{indexMatrix} but also
  \code{\linkS4class{nMatrix}}, the virtual class of boolean matrices,
  and \code{\linkS4class{sparseMatrix}}, the virtual class of sparse
  format matrices.  It is extended by \code{pMatrix}, which represents
  permutation matrices.  Objects of class \code{pMatrix} have no
  additional slots.  However, to be valid, their \code{perm} slot must
  be a permutation of \code{seq_len(n)} for \code{n=Dim[1]=Dim[2]}.
}
\seealso{
  Related classes
  \code{\linkS4class{generalMatrix}},
  \code{\linkS4class{symmetricMatrix}},
  \code{\linkS4class{posdefMatrix}},
  \code{\linkS4class{triangularMatrix}}, and
  \code{\linkS4class{diagonalMatrix}} representing
  general (unstructured, possibly nonsquare), Hermitian or symmetric,
  positive semidefinite, triangular, and diagonal matrices.
  Function \code{\link{isPerm}} for testing if an integer vector is a
  valid permutation vector.  Function \code{\link{invertPerm}} for
  computing the inverse of a permutation vector.
  \code{\link{matmult-methods}} for matrix products involving
  \code{indexMatrix}.
}
\examples{
%% moved from ./indMatrix-class.Rd :
p1 <- as(c(2,3,1), "pMatrix")
(sm1 <- as(rep(c(2,3,1), e=3), "indMatrix"))
stopifnot(all(sm1 == p1[rep(1:3, each=3),]))

## row-indexing of a <pMatrix> turns it into an <indMatrix>:
class(p1[rep(1:3, each=3),])

set.seed(12) # so we know '10' is in sample
## random index matrix for 30 observations and 10 unique values:
(s10 <- as(sample(10, 30, replace=TRUE),"indMatrix"))

## Sample rows of a numeric matrix :
(mm <- matrix(1:10, nrow=10, ncol=3))
s10 \%*\% mm

set.seed(27)
IM1 <- as(sample(1:20, 100, replace=TRUE), "indMatrix")
IM2 <- as(sample(1:18, 100, replace=TRUE), "indMatrix")
(c12 <- crossprod(IM1,IM2))
## same as cross-tabulation of the two index vectors:
stopifnot(all(c12 - unclass(table(IM1@perm, IM2@perm)) == 0))

# 3 observations, 4 implied values, first does not occur in sample:
as(2:4, "indMatrix")
# 3 observations, 5 values, first and last do not occur in sample:
as(list(2:4, 5), "indMatrix")

as(sm1, "nMatrix")
s10[1:7, 1:4] # gives an "ngTMatrix" (most economic!)
s10[1:4, ]  # preserves "indMatrix"-class

I1 <- as(c(5:1,6:4,7:3), "indMatrix")
I2 <- as(7:1, "pMatrix")
(I12 <- rbind(I1, I2))
stopifnot(is(I12, "indMatrix"),
          identical(I12, rbind(I1, I2)),
          colSums(I12) == c(2L,2:4,4:2))

%% moved from ./pMatrix-class.Rd :
(pm1 <- as(as.integer(c(2,3,1)), "pMatrix"))
t(pm1) # is the same as
solve(pm1)
pm1 \%*\% t(pm1) # check that the transpose is the inverse
stopifnot(all(diag(3) == as(pm1 \%*\% t(pm1), "matrix")),
          is.logical(as(pm1, "matrix")))

set.seed(11)
## random permutation matrix :
(p10 <- as(sample(10),"pMatrix"))

## Permute rows / columns of a numeric matrix :
(mm <- round(array(rnorm(3 * 3), c(3, 3)), 2))
mm \%*\% pm1
pm1 \%*\% mm
try(as(as.integer(c(3,3,1)), "pMatrix"))# Error: not a permutation

as(pm1, "TsparseMatrix")
p10[1:7, 1:4] # gives an "ngTMatrix" (most economic!)

## row-indexing of a <pMatrix> keeps it as an <indMatrix>:
p10[1:3, ]
}
