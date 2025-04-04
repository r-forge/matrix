\name{bandSparse}
\title{Construct Sparse Banded Matrix from (Sup-/Super-) Diagonals}
%
\keyword{array}
\keyword{utilities}
%
\alias{bandSparse}
%
\description{
  Construct a sparse banded matrix by specifying its non-zero sup- and
  super-diagonals.
}
\usage{
bandSparse(n, m = n, k, diagonals, symmetric = FALSE,
           repr = "C", giveCsparse = (repr == "C"), int2dbl = )
}
\arguments{
  \item{n,m}{the matrix dimension \eqn{(n,m) = (nrow, ncol)}.}
  \item{k}{integer vector of \dQuote{diagonal numbers}, with identical
    meaning as in \code{\link{band}(*, k)}, i.e., relative to the main diagonal,
    which is \code{k=0}.}
  \item{diagonals}{optional list of sub-/super- diagonals; if missing,
    the result will be a patter\bold{n} matrix, i.e., inheriting from
    class \code{\linkS4class{nMatrix}}.

    \code{diagonals} can also be \eqn{n' \times d}{n' x d} matrix, where
    \code{d <- length(k)} and \eqn{n' >= min(n,m)}.  In that case,
    the sub-/super- diagonals are taken from the columns of
    \code{diagonals}, where only the first several rows will be used
    (typically) for off-diagonals.
  }
  \item{symmetric}{logical; if true the result will be symmetric
    (inheriting from class \code{\linkS4class{symmetricMatrix}}) and
    only the upper or lower triangle must be specified (via \code{k} and
    \code{diagonals}).}
  \item{repr}{\code{\link{character}} string, one of \code{"C"},
    \code{"T"}, or \code{"R"}, specifying the sparse \emph{repr}esentation to
    be used for the result, i.e., one from the super classes
    \code{\linkS4class{CsparseMatrix}}, \code{\linkS4class{TsparseMatrix}}, or
    \code{\linkS4class{RsparseMatrix}}.}
  \item{giveCsparse}{(\bold{deprecated}, replaced with \code{repr}):
    logical indicating if the result should be a
    \code{\linkS4class{CsparseMatrix}} or a
    \code{\linkS4class{TsparseMatrix}}, where the default was \code{TRUE},
    and now is determined from \code{repr}; very often Csparse matrices are
    more efficient subsequently, but not always.}
  \item{int2dbl}{\code{\link{logical}} indicating if integer matrices
    (\code{\linkS4class{iMatrix}}) should be coerced to double
    (\code{\linkS4class{dMatrix}}) for back compatibility.  It is planned
    that the this will become \code{FALSE}, at first by default;
    see \code{\link{Matrix}}'s \code{int2dbl} for more details.}
}
% \details{  __needed ?__
%
% }
\value{
  a sparse matrix (of \code{\link{class}}
  \code{\linkS4class{CsparseMatrix}}) of dimension \eqn{n \times m}{n x m}
  with diagonal \dQuote{bands} as specified.
}
%\author{Martin Maechler}
\seealso{
  \code{\link{band}}, for \emph{extraction} of matrix bands;
  \code{\link{bdiag}}, \code{\link{diag}},
  \code{\link{sparseMatrix}},
  \code{\link{Matrix}}.
}
\examples{
\dontshow{ % for R_DEFAULT_PACKAGES=NULL
library(stats, pos = "package:base", verbose = FALSE)
}
diags <- list(1:30, 10*(1:20), 100*(1:20))
s1 <- bandSparse(13, k = -c(0:2, 6), diag = c(diags, diags[2]), symm=TRUE)
s1
s2 <- bandSparse(13, k =  c(0:2, 6), diag = c(diags, diags[2]), symm=TRUE)
stopifnot(identical(s1, t(s2)), is(s1,"dsCMatrix"))

## a pattern Matrix of *full* (sub-)diagonals:
bk <- c(0:4, 7,9)
(s3 <- bandSparse(30, k = bk, symm = TRUE))

## If you want a pattern matrix, but with "sparse"-diagonals,
## you currently need to go via logical sparse:
lLis <- lapply(list(rpois(20, 2), rpois(20, 1), rpois(20, 3))[c(1:3, 2:3, 3:2)],
               as.logical)
(s4 <- bandSparse(20, k = bk, symm = TRUE, diag = lLis))
(s4. <- as(drop0(s4), "nsparseMatrix"))

n <- 1e4
bk <- c(0:5, 7,11)
bMat <- matrix(1:8, n, 8, byrow=TRUE)
bLis <- as.data.frame(bMat)
B  <- bandSparse(n, k = bk, diag = bLis)
Bs <- bandSparse(n, k = bk, diag = bLis, symmetric=TRUE)
B [1:15, 1:30]
Bs[1:15, 1:30]
## can use a list *or* a matrix for specifying the diagonals:
stopifnot(identical(B,  bandSparse(n, k = bk, diag = bMat)),
	  identical(Bs, bandSparse(n, k = bk, diag = bMat, symmetric=TRUE))
          , inherits(B, "dtCMatrix") # triangular!
)
}
